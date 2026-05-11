import pandas as pd
import numpy as np
from scipy.signal import medfilt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- Parameters ---
FILE_COUNT = 30
VELOCITY_THRESHOLD = 0.0085
COALESCENCE_TIME = 0.135
RINGING_FREQ_HZ = 6

# --- Threshold candidates to overlay ---
FIXED_THRESHOLDS = [1.0, 2.0, 5.0 ,10.0]  # seconds
THRESHOLD_COLORS = ["blue", "green", "orange", "red"]

# --- Build subplots: 10 rows x 3 cols ---
n_cols = 3
n_rows = 10

fig = make_subplots(
    rows=n_rows, cols=n_cols,
    subplot_titles=[f"Run {i}" for i in range(1, FILE_COUNT + 1)],
    horizontal_spacing=0.06,
    vertical_spacing=0.04
)

for i in range(1, FILE_COUNT + 1):
    filename = f"{i}.csv"
    row = (i - 1) // n_cols + 1
    col = (i - 1) % n_cols + 1
    show_legend = (i == 1)

    try:
        df = pd.read_csv(
            filename,
            names=["t_s", "count", "theta_rad", "x_m"],
            comment="#",
            on_bad_lines="skip",
            dtype=str
        )
        for c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.dropna()

        dt_median = df["t_s"].diff().median()
        if pd.isna(dt_median) or dt_median == 0:
            continue
        fs = 1 / dt_median
        window_samples = int((1 / RINGING_FREQ_HZ) * fs)
        if window_samples % 2 == 0:
            window_samples += 1

        df["x_filtered"] = medfilt(df["x_m"], kernel_size=window_samples)
        df["v_filtered"] = df["x_filtered"].diff() / df["t_s"].diff()
        df["v_filtered_abs"] = df["v_filtered"].abs()

    except Exception as e:
        print(f"Skipping run {i}: {e}")
        continue

    # Load events for this run
    events_df = pd.read_csv("event_data_shifted.csv")
    events = events_df[events_df["file_id"] == i].sort_values("time").reset_index(drop=True)

    if len(events) < 2:
        continue

    # Compute waiting times from event times
    times = events["time"].values
    wait_times = np.diff(times)
    wait_t = times[1:]  # time of the second event in each pair

    # Plot Δt vs t (the aging plot per run)
    fig.add_trace(go.Scatter(
        x=wait_t,
        y=wait_times,
        mode="markers",
        marker=dict(size=4, color="gray", opacity=0.7),
        name="Δt",
        showlegend=False
    ), row=row, col=col)

    # Overlay fixed threshold lines
    t_range = np.array([times.min(), times.max()])
    for thresh, color in zip(FIXED_THRESHOLDS, THRESHOLD_COLORS):
        fig.add_trace(go.Scatter(
            x=t_range,
            y=[thresh, thresh],
            mode="lines",
            line=dict(color=color, width=1.5, dash="dot"),
            name=f"Δt={thresh}s",
            showlegend=show_legend,
            legendgroup=f"thresh_{thresh}"
        ), row=row, col=col)
        show_legend = False  # only show legend once per threshold

    fig.update_xaxes(type="log", row=row, col=col)
    fig.update_yaxes(type="log", row=row, col=col)

# Add legend manually for thresholds
for thresh, color in zip(FIXED_THRESHOLDS, THRESHOLD_COLORS):
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode="lines",
        line=dict(color=color, width=2, dash="dot"),
        name=f"fixed Δt = {thresh}s",
        showlegend=True
    ))

fig.update_layout(
    title_text="Per-run Δt vs t — fixed threshold candidates (log-log)",
    height=2500,
    width=1400,
    template="plotly_white",
    legend=dict(
        x=1.01, y=1,
        bordercolor="black",
        borderwidth=1
    )
)

fig.show()
# fig.write_html("threshold_verification_all_runs.html")