import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

CSV = "7.2HZ.csv"
VELOCITY_THRESHOLD = 0.005  # <--- ADJUST THIS BASED ON THE HISTOGRAM

# --- Data Loading & Cleaning ---
df = pd.read_csv(
    CSV,
    names=["t_s", "count", "theta_rad", "x_m"],
    comment="#",
    on_bad_lines="skip"
)

# Convert to numeric and drop errors
for c in df.columns:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna()

# --- Calculate Velocity ---
# We use .abs() for the threshold check, but plot the real sign
df["v_mps"] = df["x_m"].diff() / df["t_s"].diff()
df["v_abs"] = df["v_mps"].abs()

# Identify "Events" (Where velocity > threshold)
df["is_event"] = df["v_abs"] > VELOCITY_THRESHOLD

# --- Plotting ---
fig = make_subplots(
    rows=2, cols=2,
    specs=[[{"colspan": 2}, None], [{"colspan": 1}, {"colspan": 1}]],
    column_widths=[0.7, 0.3],
    subplot_titles=("Velocity Time Series (Red = Event)", "Displacement", "Velocity Distribution"),
    vertical_spacing=0.15
)

# 1. Velocity vs Time (Color-coded by threshold)
# Background (Noise)
fig.add_trace(go.Scattergl(
    x=df.loc[~df["is_event"], "t_s"],
    y=df.loc[~df["is_event"], "v_mps"],
    mode="markers",
    marker=dict(color='gray', size=2, opacity=0.5),
    name="Noise (< Threshold)"
), row=1, col=1)

# Foreground (Events)
fig.add_trace(go.Scattergl(
    x=df.loc[df["is_event"], "t_s"],
    y=df.loc[df["is_event"], "v_mps"],
    mode="markers",
    marker=dict(color='red', size=4),
    name="Slip Event (> Threshold)"
), row=1, col=1)

# Add the Threshold Line
fig.add_hline(y=VELOCITY_THRESHOLD, line_dash="dash", line_color="black", row=1, col=1, annotation_text="Threshold")
fig.add_hline(y=-VELOCITY_THRESHOLD, line_dash="dash", line_color="black", row=1, col=1)

# 2. Displacement vs Time (Reference)
fig.add_trace(go.Scattergl(
    x=df["t_s"], y=df["x_m"],
    mode="lines", line=dict(color='blue', width=1),
    name="Displacement"
), row=2, col=1)

# 3. Histogram of Log Velocity (The Decision Maker)
# We filter out zeros to avoid log(0) errors
v_non_zero = df[df["v_abs"] > 1e-6]["v_abs"]
fig.add_trace(go.Histogram(
    x=np.log10(v_non_zero),
    nbinsx=50,
    marker_color='black',
    name="Log10(|v|) Hist"
), row=2, col=2)

# Add Threshold Line to Histogram
fig.add_vline(x=np.log10(VELOCITY_THRESHOLD), line_dash="dash", line_color="red", row=2, col=2, annotation_text="Cutoff")

# Update Layout
fig.update_layout(
    height=800,
    title_text=f"Event Detection (Threshold = {VELOCITY_THRESHOLD} m/s)",
    hovermode="closest",
    showlegend=True
)

fig.update_xaxes(title_text="Time [s]", row=1, col=1)
fig.update_xaxes(title_text="Time [s]", row=2, col=1)
fig.update_xaxes(title_text="Log10(Velocity)", row=2, col=2)
fig.update_yaxes(title_text="Velocity [m/s]", row=1, col=1)
fig.update_yaxes(title_text="Displacement [m]", row=2, col=1)
fig.update_yaxes(title_text="Count", row=2, col=2)

fig.show()
