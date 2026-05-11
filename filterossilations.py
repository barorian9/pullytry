import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from scipy.signal import medfilt

# --- Parameters ---
CSV = "24.csv"
VELOCITY_THRESHOLD = 0.0085   # Threshold for the SMOOTHED velocity
COALESCENCE_TIME = 0.065     # Merge events closer than this (seconds)
RINGING_FREQ_HZ = 6       # The frequency to filter out

# --- Load & Prep ---
df = pd.read_csv(CSV, names=["t_s", "count", "theta_rad", "x_m"], comment="#", on_bad_lines="skip")
for c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna()

# --- 1. FILTERING (The Crucial Step) ---
# Calculate Sampling Rate
dt_median = df["t_s"].diff().median()
fs = 1 / dt_median

# Calculate Window Size (1 Period of the Ringing)
period_sec = 1 / RINGING_FREQ_HZ
window_samples = int(period_sec * fs)
if window_samples % 2 == 0: window_samples += 1  # Must be odd

print(f"Filter: Removing {RINGING_FREQ_HZ}Hz (Window: {window_samples} samples)")

# Apply Median Filter to DISPLACEMENT
df["x_filtered"] = medfilt(df["x_m"], kernel_size=window_samples)

# Calculate Velocity from the FILTERED Displacement
# (This effectively 'smooths' the velocity too)
df["v_filtered"] = df["x_filtered"].diff() / df["t_s"].diff()
df["v_filtered_abs"] = df["v_filtered"].abs()

# --- 2. EVENT DETECTION (On Smoothed Data) ---
# We use 'v_filtered_abs' here, so we aren't detecting the ringing noise
df["is_active"] = df["v_filtered"] < -VELOCITY_THRESHOLD

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

print(f"Detected {len(events_merged)} events on smoothed data.")

# --- 3. PLOTTING ---
fig = make_subplots(
    rows=2, cols=1,
    shared_xaxes=True,
    vertical_spacing=0.05,
    subplot_titles=(f"Displacement (Raw vs Filtered @ {RINGING_FREQ_HZ}Hz)", "Velocity (Derived from Filtered Data)")
)

# --- ROW 1: Displacement ---
# Raw (Background, Faint)
fig.add_trace(go.Scatter(
    x=df["t_s"], y=df["x_m"],
    mode="lines",
    line=dict(color='red', width=1),
    opacity=0.5,
    name="Raw Displacement"
), row=1, col=1)

# Filtered (Foreground, Dark)
fig.add_trace(go.Scatter(
    x=df["t_s"], y=df["x_filtered"],
    mode="lines",
    line=dict(color='#2c3e50', width=2),
    name="Filtered Displacement"
), row=1, col=1)

# --- ROW 2: Velocity ---
# This is the velocity of the CLEAN signal
fig.add_trace(go.Scatter(
    x=df["t_s"], y=df["v_filtered"],
    mode="lines",
    line=dict(color='#e74c3c', width=1),
    name="Filtered Velocity"
), row=2, col=1)

# Plot the Threshold Crossings (on the filtered data)
fig.add_trace(go.Scatter(
    x=df.loc[df["is_active"], "t_s"],
    y=df.loc[df["is_active"], "v_filtered"],
    mode="markers",
    marker=dict(color='black', size=4),
    name="Active (Filtered)"
), row=2, col=1)

# --- Event Highlights ---
for start, end in events_merged:
    if start == end: end += 0.05
    fig.add_vrect(
        x0=start, x1=end,
        fillcolor="green", opacity=0.1,
        layer="below", line_width=0
    )

fig.update_layout(
    title=f"Analysis of Smoothed Data (Filter: {RINGING_FREQ_HZ}Hz)",
    hovermode="x unified",
    height=800,
    showlegend=True
)

fig.update_yaxes(title_text="Position [m]", row=1, col=1)
fig.update_yaxes(title_text="Velocity [m/s]", row=2, col=1)
fig.update_xaxes(title_text="Time [s]", row=2, col=1)

fig.show()