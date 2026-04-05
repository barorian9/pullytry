import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt

# --- Parameters ---
CSV = "7.2HZ.csv"
VELOCITY_THRESHOLD = 0.0085  # Threshold for filtered velocity
COALESCENCE_TIME = 0.135  # Merge events closer than this (s)
RINGING_FREQ_HZ = 6  # Ringing frequency to filter out

# --- 1. Load & Filter Data ---
df = pd.read_csv(CSV, names=["t_s", "count", "theta_rad", "x_m"], comment="#", on_bad_lines="skip")
for c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna()

# Calculate Filter Window (1 Period of Ringing)
dt_median = df["t_s"].diff().median()
fs = 1 / dt_median
period_sec = 1 / RINGING_FREQ_HZ
window_samples = int(period_sec * fs)
if window_samples % 2 == 0: window_samples += 1

print(f"Applying Median Filter (Window: {window_samples} samples)")
df["x_filtered"] = medfilt(df["x_m"], kernel_size=window_samples)

# Calculate Velocity from SMOOTHED displacement
df["v_filtered"] = df["x_filtered"].diff() / df["t_s"].diff()
df["v_filtered_abs"] = df["v_filtered"].abs()

# --- 2. Event Detection ---
df["is_active"] = df["v_filtered_abs"] > VELOCITY_THRESHOLD
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

print(f"Detected {len(events_merged)} events.")

# --- 3. Calculate Metrics (t, dt, Magnitude) ---
event_metrics = []

for i, (start, end) in enumerate(events_merged):
    # Find indices for this event window
    mask = (df["t_s"] >= start) & (df["t_s"] <= end)
    event_data = df[mask]

    if event_data.empty: continue

    # A. Magnitude: Displacement from Start to End
    # We use x_filtered to avoid endpoint noise
    x_start = df.loc[df["t_s"] >= start, "x_filtered"].iloc[0]
    x_end = df.loc[df["t_s"] <= end, "x_filtered"].iloc[-1]
    magnitude = abs(x_end - x_start)

    # B. Exact Time (t_peak): Time of max velocity
    t_peak = event_data.loc[event_data["v_filtered_abs"].idxmax(), "t_s"]

    # C. Waiting Time (dt): Time since PREVIOUS event ended
    # (Start of current - End of previous)
    if i == 0:
        waiting_time = np.nan
    else:
        prev_end = events_merged[i - 1][1]
        waiting_time = start - prev_end

    event_metrics.append({
        "time": t_peak,
        "waiting_time": waiting_time,
        "magnitude": magnitude,
        "duration": end - start
    })

# Save to DataFrame
results_df = pd.DataFrame(event_metrics)

# --- 4. Plotting ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Plot 1: Waiting Time vs. Time (The "Aging" Plot)
# Log-Log is standard to check for power-law growth
axes[0].semilogy(results_df["time"], results_df["waiting_time"], 'o',
                 color='#2c3e50', alpha=0.6, markersize=4)
axes[0].set_xlabel("Time $t$ (s) [System Age]", fontsize=12)
axes[0].set_ylabel("Waiting Time $\Delta t$ (s)", fontsize=12)
axes[0].set_title("Aging: $\Delta t$ vs $t$", fontsize=14)
axes[0].grid(True, which="both", alpha=0.3)

# Plot 2: Histogram of Waiting Times
axes[1].hist(results_df["waiting_time"].dropna(), bins=30,
             color='#16a085', edgecolor='black', alpha=0.7)
axes[1].set_xlabel("Waiting Time $\Delta t$ (s)", fontsize=12)
axes[1].set_ylabel("Count", fontsize=12)
axes[1].set_title("Distribution of Waiting Times", fontsize=14)
axes[1].set_yscale('log')  # Log scale y-axis often reveals rare long waits
axes[1].grid(True, alpha=0.3)

# Plot 3: Histogram of Magnitudes (Gutenberg-Richter check)
axes[2].hist(results_df["magnitude"], bins=30,
             color='#8e44ad', edgecolor='black', alpha=0.7)
axes[2].set_xlabel("Magnitude $\Delta x$ (m)", fontsize=12)
axes[2].set_ylabel("Count", fontsize=12)
axes[2].set_title("Distribution of Magnitudes", fontsize=14)
axes[2].set_yscale('log')  # Gutenberg-Richter implies a straight line here
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Optional: Print summary
print("Top 5 Largest Avalanches:")
print(results_df.nlargest(5, "magnitude")[["time", "magnitude", "waiting_time"]])