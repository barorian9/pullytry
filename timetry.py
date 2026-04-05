import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

# --- Parameters ---
CSV = "7.2HZ.csv"
VELOCITY_THRESHOLD = 0.012
COALESCENCE_TIME = 0.135  # <--- TWEAK THIS (Seconds). Try 0.1, 0.5, 2.0

# --- Load & Prep ---
df = pd.read_csv(CSV, names=["t_s", "count", "theta_rad", "x_m"], comment="#", on_bad_lines="skip")
for c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna()

# Calculate Velocity
df["v_mps"] = df["x_m"].diff() / df["t_s"].diff()
df["v_abs"] = df["v_mps"].abs()

# 1. Identify Raw Threshold Crossings
df["is_active"] = df["v_abs"] > VELOCITY_THRESHOLD

# 2. Merging Logic
active_indices = df[df["is_active"]].index
events_merged = []

if len(active_indices) > 0:
    active_times = df.loc[active_indices, "t_s"].values
    time_diffs = np.diff(active_times)
    break_points = np.where(time_diffs > COALESCENCE_TIME)[0]

    event_starts = [active_times[0]]
    event_ends = []

    for bp in break_points:
        event_ends.append(active_times[bp])
        event_starts.append(active_times[bp + 1])

    event_ends.append(active_times[-1])
    events_merged = list(zip(event_starts, event_ends))

# --- Plotting ---
# Create a subplot with 2 rows (Displacement top, Velocity bottom)
fig = make_subplots(
    rows=2, cols=1,
    shared_xaxes=True,
    vertical_spacing=0.05,
    subplot_titles=("Displacement (x vs t)", "Velocity (v vs t)")
)

# --- ROW 1: Displacement (x vs t) ---
fig.add_trace(go.Scattergl(
    x=df["t_s"], y=df["x_m"],
    mode="lines",
    line=dict(color='#2c3e50', width=2),
    name="Displacement"
), row=1, col=1)

# --- ROW 2: Velocity (v vs t) ---
# Raw Velocity Line
fig.add_trace(go.Scattergl(
    x=df["t_s"], y=df["v_mps"],
    mode="lines",
    line=dict(color='lightgray', width=1),
    name="Velocity"
), row=2, col=1)

# Threshold Crossings (Red Dots)
fig.add_trace(go.Scattergl(
    x=df.loc[df["is_active"], "t_s"],
    y=df.loc[df["is_active"], "v_mps"],
    mode="markers",
    marker=dict(color='red', size=4),
    name="Raw Crossings"
), row=2, col=1)

# --- Event Highlights (Vertical Rectangles) ---
# This adds the green shaded region to BOTH subplots automatically
for start, end in events_merged:
    if start == end: end += 0.05

    fig.add_vrect(
        x0=start, x1=end,
        fillcolor="green", opacity=0.1,
        layer="below", line_width=0
    )

# --- Layout Polish ---
fig.update_layout(
    title=f"Event Analysis (Wait Time: {COALESCENCE_TIME}s)",
    hovermode="x unified",
    height=800,  # Taller figure to fit both graphs
    showlegend=True
)

# Axis Labels
fig.update_yaxes(title_text="Position [m]", row=1, col=1)
fig.update_yaxes(title_text="Velocity [m/s]", row=2, col=1)
fig.update_xaxes(title_text="Time [s]", row=2, col=1)

fig.show()