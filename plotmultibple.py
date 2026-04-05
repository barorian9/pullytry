import pandas as pd
import numpy as np
from scipy.signal import medfilt
import os
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- Parameters ---
FILE_COUNT = 30
VELOCITY_THRESHOLD = 0.0085
COALESCENCE_TIME = 0.135
RINGING_FREQ_HZ = 6

all_event_metrics = []

# --- 1. Data Processing Loop (Unchanged) ---
for i in range(1, FILE_COUNT + 1):
    filename = f"{i}.csv"
    if not os.path.exists(filename):
        print(f"Skipping {filename}: File not found.")
        continue

    print(f"Processing {filename}...")

    df = pd.read_csv(
        filename,
        names=["t_s", "count", "theta_rad", "x_m"],
        comment="#",
        on_bad_lines="skip",
        dtype=str
    )
    for c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna()

    dt_median = df["t_s"].diff().median()
    if pd.isna(dt_median) or dt_median == 0: continue  # Safety check

    fs = 1 / dt_median
    window_samples = int((1 / RINGING_FREQ_HZ) * fs)
    if window_samples % 2 == 0: window_samples += 1

    df["x_filtered"] = medfilt(df["x_m"], kernel_size=window_samples)
    df["v_filtered"] = df["x_filtered"].diff() / df["t_s"].diff()
    df["v_filtered_abs"] = df["v_filtered"].abs()

    df["is_active"] = df["v_filtered_abs"] > VELOCITY_THRESHOLD
    active_indices = df[df["is_active"]].index

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

        for j, (start, end) in enumerate(events_merged):
            mask = (df["t_s"] >= start) & (df["t_s"] <= end)
            event_data = df[mask]
            if event_data.empty: continue

            x_start = df.loc[df["t_s"] >= start, "x_filtered"].iloc[0]
            x_end = df.loc[df["t_s"] <= end, "x_filtered"].iloc[-1]

            all_event_metrics.append({
                "file_id": i,
                "time": event_data.loc[event_data["v_filtered_abs"].idxmax(), "t_s"],
                "waiting_time": start - events_merged[j - 1][1] if j > 0 else np.nan,
                "magnitude": abs(x_end - x_start)
            })

# --- 2. Interactive Plotting with Plotly ---
results_df = pd.DataFrame(all_event_metrics)

if results_df.empty:
    print("No events detected. Check threshold parameters.")
else:
    # Create a subplot figure with 1 row and 3 columns
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=("Combined Aging Plot", "Global Waiting Times", "Global Magnitudes")
    )

    # Plot 1: Aging (Scatter Plot)
    # Using Scattergl for performance with large datasets
    fig.add_trace(
        go.Scattergl(
            x=results_df["time"],
            y=results_df["waiting_time"],
            mode='markers',
            marker=dict(
                size=5,
                color=results_df["file_id"],  # Color by file ID
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(title="File ID", x=-0.1)  # Move colorbar to left
            ),
            text=results_df["file_id"],  # Shows file ID on hover
            hovertemplate="<b>File: %{text}</b><br>Time: %{x:.2f}s<br>Wait: %{y:.4f}s<extra></extra>",
            name="Events"
        ),
        row=1, col=1
    )

    # Plot 2: Waiting Time Distribution (Histogram)
    fig.add_trace(
        go.Histogram(
            x=results_df["waiting_time"].dropna(),
            nbinsx=1000,
            marker_color='#16a085',
            name="Wait Dist",
            hovertemplate="Wait: %{x:.4f}s<br>Count: %{y}<extra></extra>"
        ),
        row=1, col=2
    )

    # Plot 3: Magnitude Distribution (Histogram)
    fig.add_trace(
        go.Histogram(
            x=results_df["magnitude"],
            nbinsx=1000,
            marker_color='#8e44ad',
            name="Mag Dist",
            hovertemplate="Mag: %{x:.6f}m<br>Count: %{y}<extra></extra>"
        ),
        row=1, col=3
    )

    # --- Layout Updates ---

    # Update Axes for Plot 1 (Aging)
    fig.update_xaxes(title_text="Time t (s)", row=1, col=1)
    fig.update_yaxes(title_text="Waiting Time Δt (s)", type="log", row=1, col=1)

    # Update Axes for Plot 2 (Wait Hist)
    fig.update_xaxes(title_text="Waiting Time (s)", row=1, col=2)
    fig.update_yaxes(title_text="Count", type="log", row=1, col=2)

    # Update Axes for Plot 3 (Mag Hist)
    fig.update_xaxes(title_text="Magnitude (m)", row=1, col=3)
    fig.update_yaxes(title_text="Count", type="log", row=1, col=3)

    # General Layout
    fig.update_layout(
        title_text="Interactive Event Analysis",
        showlegend=False,
        height=500,
        width=1500,
        template="plotly_white"  # Cleaner look
    )

    # Show and Save
    fig.show()
    # fig.write_html("interactive_results.html") # Uncomment to save as file