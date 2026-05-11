import pandas as pd
import numpy as np
from scipy.signal import medfilt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- Parameters (must match the code that produced event_data.csv) ---
RUN = 1                    # Which run to inspect — change this
CSV = f"{RUN}.csv"
VELOCITY_THRESHOLD = 0.0085
COALESCENCE_TIME = 0.135
RINGING_FREQ_HZ = 6

# Optional: zoom into a time window (set to None to see full run)
T_MIN = None   # e.g. 50
T_MAX = None   # e.g. 200

# --- Load raw data ---
df = pd.read_csv(
    CSV,
    names=["t_s", "count", "theta_rad", "x_m"],
    comment="#",
    on_bad_lines="skip",
    dtype=str
)
for c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna()

# --- Filter ---
dt_median = df["t_s"].diff().median()
fs = 1 / dt_median
window_samples = int((1 / RINGING_FREQ_HZ) * fs)
if window_samples % 2 == 0: window_samples += 1

df["x_filtered"] = medfilt(df["x_m"], kernel_size=window_samples)
df["v_filtered"] = df["x_filtered"].diff() / df["t_s"].diff()
df["v_filtered_abs"] = df["v_filtered"].abs()

# --- Load detected events from event_data.csv ---
events = pd.read_csv("event_data.csv")
events = events[events["file_id"] == RUN].sort_values("time").reset_index(drop=True)

print(f"Run {RUN}: {len(df)} raw samples, {len(events)} detected events")
print(f"Sampling rate: {fs:.1f} Hz")
print(f"Filter window: {window_samples} samples ({1/RINGING_FREQ_HZ:.3f}s)")

# --- Apply zoom if set ---
if T_MIN is not None:
    df = df[df["t_s"] >= T_MIN]
    events = events[events["time"] >= T_MIN]
if T_MAX is not None:
    df = df[df["t_s"] <= T_MAX]
    events = events[events["time"] <= T_MAX]

# --- Plot ---
fig = make_subplots(
    rows=2, cols=1,
    shared_xaxes=True,
    vertical_spacing=0.05,
    subplot_titles=(
        f"Run {RUN} — Displacement (zoom T={T_MIN}-{T_MAX}s)",
        "Filtered Velocity"
    )
)

# Row 1: Raw displacement (faint)
fig.add_trace(go.Scatter(
    x=df["t_s"], y=df["x_m"],
    mode="lines",
    line=dict(color="red", width=0.8),
    opacity=0.4,
    name="Raw displacement"
), row=1, col=1)

# Row 1: Filtered displacement
fig.add_trace(go.Scatter(
    x=df["t_s"], y=df["x_filtered"],
    mode="lines",
    line=dict(color="black", width=1.5),
    name="Filtered displacement"
), row=1, col=1)

# Row 1: Event markers on displacement
fig.add_trace(go.Scatter(
    x=events["time"],
    y=np.interp(events["time"], df["t_s"].values, df["x_filtered"].values),
    mode="markers",
    marker=dict(color="green", size=8, symbol="triangle-down"),
    name="Detected events"
), row=1, col=1)

# Row 2: Filtered velocity
fig.add_trace(go.Scatter(
    x=df["t_s"], y=df["v_filtered"],
    mode="lines",
    line=dict(color="#e74c3c", width=1),
    name="Filtered velocity"
), row=2, col=1)

# Row 2: Threshold lines
fig.add_hline(
    y=VELOCITY_THRESHOLD,
    line=dict(color="blue", width=1, dash="dash"),
    row=2, col=1
)
fig.add_hline(
    y=-VELOCITY_THRESHOLD,
    line=dict(color="blue", width=1, dash="dash"),
    row=2, col=1
)

# Row 2: Active samples
active = df[df["v_filtered_abs"] > VELOCITY_THRESHOLD]
fig.add_trace(go.Scatter(
    x=active["t_s"],
    y=active["v_filtered"],
    mode="markers",
    marker=dict(color="black", size=3),
    name="Above threshold"
), row=2, col=1)

# Vertical lines at each detected event
for _, ev in events.iterrows():
    fig.add_vline(
        x=ev["time"],
        line=dict(color="green", width=1, dash="dot"),
        row="all", col=1
    )

fig.update_yaxes(title_text="Position [m]", row=1, col=1)
fig.update_yaxes(title_text="Velocity [m/s]", row=2, col=1)
fig.update_xaxes(title_text="Time [s]", row=2, col=1)
fig.update_layout(
    height=700,
    width=1200,
    template="plotly_white",
    hovermode="x unified",
    title_text=f"Event detection verification — Run {RUN}"
)

fig.show()