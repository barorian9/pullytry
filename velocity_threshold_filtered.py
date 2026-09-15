import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.signal import medfilt
import numpy as np

CSV = "7.2HZ.csv"
VELOCITY_THRESHOLD = 0.0085  # pipeline value, being tested here
RINGING_FREQ_HZ = 6          # A4/A6: filter kernel = floor(fs/RINGING_FREQ_HZ), forced odd

# --- Data Loading & Cleaning (same as derivative_threshold.py) ---
df = pd.read_csv(
    CSV,
    names=["t_s", "count", "theta_rad", "x_m"],
    comment="#",
    on_bad_lines="skip"
)
for c in df.columns:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna()

# --- Median filter (A4/A6: kernel=17 at fs=100Hz) ---
dt_median = df["t_s"].diff().median()
fs = 1 / dt_median
window_samples = int((1 / RINGING_FREQ_HZ) * fs)
if window_samples % 2 == 0:
    window_samples += 1
print(f"fs = {fs:.3f} Hz, filter kernel = {window_samples} samples")

df["x_filtered"] = medfilt(df["x_m"], kernel_size=window_samples)

# --- Calculate velocity on FILTERED displacement (pipeline's actual input) ---
df["v_mps"] = df["x_filtered"].diff() / df["t_s"].diff()
df["v_abs"] = df["v_mps"].abs()
df["is_event"] = df["v_abs"] > VELOCITY_THRESHOLD

# --- Plotting (mirrors derivative_threshold.py, filtered data) ---
fig = make_subplots(
    rows=2, cols=2,
    specs=[[{"colspan": 2}, None], [{"colspan": 1}, {"colspan": 1}]],
    column_widths=[0.7, 0.3],
    subplot_titles=("Filtered Velocity Time Series (Red = Event)", "Filtered Displacement", "Filtered Velocity Distribution"),
    vertical_spacing=0.15
)

fig.add_trace(go.Scatter(
    x=df.loc[~df["is_event"], "t_s"],
    y=df.loc[~df["is_event"], "v_mps"],
    mode="markers",
    marker=dict(color='gray', size=2, opacity=0.5),
    name="Noise (< Threshold)"
), row=1, col=1)

fig.add_trace(go.Scatter(
    x=df.loc[df["is_event"], "t_s"],
    y=df.loc[df["is_event"], "v_mps"],
    mode="markers",
    marker=dict(color='red', size=4),
    name="Slip Event (> Threshold)"
), row=1, col=1)

fig.add_hline(y=VELOCITY_THRESHOLD, line_dash="dash", line_color="black", row=1, col=1, annotation_text="Threshold")
fig.add_hline(y=-VELOCITY_THRESHOLD, line_dash="dash", line_color="black", row=1, col=1)

fig.add_trace(go.Scatter(
    x=df["t_s"], y=df["x_filtered"],
    mode="lines", line=dict(color='blue', width=1),
    name="Filtered Displacement"
), row=2, col=1)

v_non_zero = df[df["v_abs"] > 1e-6]["v_abs"]
fig.add_trace(go.Histogram(
    x=np.log10(v_non_zero),
    nbinsx=50,
    marker_color='black',
    name="Log10(|v_filtered|) Hist"
), row=2, col=2)

fig.add_vline(x=np.log10(VELOCITY_THRESHOLD), line_dash="dash", line_color="red", row=2, col=2, annotation_text="Pipeline threshold")

fig.update_layout(
    height=800,
    title_text=f"Event Detection on FILTERED velocity (Threshold = {VELOCITY_THRESHOLD} m/s)",
    hovermode="closest",
    showlegend=True
)

fig.update_xaxes(title_text="Time [s]", row=1, col=1)
fig.update_xaxes(title_text="Time [s]", row=2, col=1)
fig.update_xaxes(title_text="Log10(Velocity)", row=2, col=2)
fig.update_yaxes(title_text="Velocity [m/s]", row=1, col=1)
fig.update_yaxes(title_text="Displacement [m]", row=2, col=1)
fig.update_yaxes(title_text="Count", row=2, col=2)

fig.write_html("velocity_threshold_filtered.html")
print("Saved velocity_threshold_filtered.html")

# --- Numeric check: where is the noise/event gap in the filtered histogram? ---
log_v = np.log10(v_non_zero.values)
counts, edges = np.histogram(log_v, bins=50)
centers = (edges[:-1] + edges[1:]) / 2

# look for the valley (local minimum) between the two modes, restricted to
# the plausible threshold range so we don't pick up tail noise
mask = (centers > -3.5) & (centers < -1.5)
if mask.sum() > 2:
    idx = np.where(mask)[0]
    valley_i = idx[np.argmin(counts[idx])]
    valley_v = 10 ** centers[valley_i]
    print(f"Histogram valley (candidate noise/event boundary) near v = {valley_v:.5f} m/s")
else:
    print("Not enough bins in search range to locate a valley")

print(f"Pipeline VELOCITY_THRESHOLD = {VELOCITY_THRESHOLD} m/s")
print(f"Fraction of samples flagged as events at this threshold: {df['is_event'].mean():.4%}")