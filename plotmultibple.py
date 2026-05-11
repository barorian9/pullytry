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

# --- Log Bin Limits --- SET YOUR LIMITS HERE ---
WAIT_MIN, WAIT_MAX = 1e-3, 100
MAG_MIN,  MAG_MAX  = 1e-6, 1e-1
N_BINS = 50

all_event_metrics = []

# --- 1. Data Processing Loop ---
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
    if pd.isna(dt_median) or dt_median == 0: continue

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
results_df.to_csv("event_data.csv", index=False)  # raw times preserved

# Shifted version: t=0 is first detected event per run
shifted_df = results_df.copy()
t0_per_run = shifted_df.groupby("file_id")["time"].transform("min")
shifted_df["time"] = shifted_df["time"] - t0_per_run
shifted_df.to_csv("event_data_shifted.csv", index=False)

if results_df.empty:
    print("No events detected. Check threshold parameters.")
else:
    aging_df = shifted_df[shifted_df["time"] > np.exp(0)]

    wait_bins = np.geomspace(WAIT_MIN, WAIT_MAX, N_BINS + 1)
    mag_bins  = np.geomspace(MAG_MIN,  MAG_MAX,  N_BINS + 1)

    wait_data = results_df["waiting_time"].dropna()
    mag_data  = results_df["magnitude"]

    wait_counts, _ = np.histogram(wait_data, bins=wait_bins, density=True)
    mag_counts,  _ = np.histogram(mag_data,  bins=mag_bins,  density=True)

    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=("Combined Aging Plot", "Global Waiting Times PDF", "Global Magnitudes PDF")
    )

    # Plot 1: Aging
    fig.add_trace(
        go.Scatter(
            x=aging_df["time"],  # raw time, no np.log
            y=aging_df["waiting_time"],
            mode='markers',
            marker=dict(
                size=5,
                color=aging_df["file_id"],
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(title="File ID", x=-0.1)
            ),
            text=aging_df["file_id"],
            hovertemplate="<b>File: %{text}</b><br>Time: %{x:.2f}s<br>Wait: %{y:.4f}s<extra></extra>",
            name="Events",
            showlegend=False
        ),
        row=1, col=1
    )




    # Plot 2: Waiting Time PDF
    fig.add_trace(
        go.Bar(
            x=wait_bins[:-1],
            y=wait_counts,
            width=np.diff(wait_bins),
            marker_color='#16a085',
            name="Wait PDF",
            showlegend=False,
            hovertemplate="Wait: %{x:.4f}s<br>PDF: %{y:.4e}<extra></extra>"
        ),
        row=1, col=2
    )

    # Omori fit: t^-1 from 0.2 to 20
    x_omori = np.geomspace(0.2, 30, 200)
    x_ref = 2.0
    y_ref = wait_counts[np.argmin(np.abs(wait_bins[:-1] - x_ref))]
    scale_omori = y_ref / (x_ref ** -1)
    y_omori = scale_omori * x_omori ** (-1)

    fig.add_trace(
        go.Scatter(
            x=x_omori,
            y=y_omori,
            mode='lines',
            line=dict(color='red', width=4, dash='dash'),
            showlegend=False,
            hovertemplate="t: %{x:.2f}s<br>power t⁻¹: %{y:.4e}<extra></extra>"
        ),
        row=1, col=2
    )

    # Power -2 tail anchored at t=20
    x_tail = np.geomspace(15, WAIT_MAX, 200)
    y_at_20 = scale_omori * (20 ** -1)
    scale_tail = y_at_20 / (20 ** -2)
    y_tail = scale_tail * x_tail ** (-2)

    fig.add_trace(
        go.Scatter(
            x=x_tail,
            y=y_tail,
            mode='lines',
            line=dict(color='orange', width=4, dash='dash'),
            showlegend=False,
            hovertemplate="t: %{x:.2f}s<br>power t⁻²: %{y:.4e}<extra></extra>"
        ),
        row=1, col=2
    )

    # Inline legend annotations inside plot 2
    for text, color, y_pos in [
        ("power ∝ t⁻¹", "red",    0.95),
        ("power ∝ t⁻²",        "orange", 0.84),
    ]:
        fig.add_annotation(
            x=0.97, y=y_pos,
            xref="x2 domain", yref="y2 domain",
            text=f"<b>— — {text}</b>",
            showarrow=False,
            font=dict(size=12, color=color),
            align="right",
            bgcolor="white",
            bordercolor=color,
            borderwidth=1,
            opacity=0.9,
        )

    # Plot 3: Magnitude PDF
    fig.add_trace(
        go.Bar(
            x=mag_bins[:-1],
            y=mag_counts,
            width=np.diff(mag_bins),
            marker_color='#8e44ad',
            name="Mag PDF",
            showlegend=False,
            hovertemplate="Mag: %{x:.6f}m<br>PDF: %{y:.4e}<extra></extra>"
        ),
        row=1, col=3
    )

    # --- Layout Updates ---
    fig.update_xaxes(title_text="Time t (s)", type="log", row=1, col=1)
    fig.update_yaxes(title_text="Waiting Time Δt (s)", type="log", row=1, col=1)

    fig.update_xaxes(title_text="Waiting Time (s)", type="log", row=1, col=2)
    fig.update_yaxes(title_text="PDF (1/s)", type="log", row=1, col=2)

    fig.update_xaxes(title_text="Magnitude (m)", type="log", row=1, col=3)
    fig.update_yaxes(title_text="PDF (1/m)", type="log", row=1, col=3)

    fig.update_layout(
        title_text="Interactive Event Analysis",
        showlegend=False,
        height=500,
        width=1600,
        template="plotly_white"
    )

    fig.show()
    # fig.write_html("interactive_results.html")