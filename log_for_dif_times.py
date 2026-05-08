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

# --- Log Bin Limits ---
WAIT_MIN, WAIT_MAX = 1e-3, 100
N_BINS = 50

# --- Time Windows ---
WINDOWS = [
    {"label": "ln(t) 2–3", "t_min": np.exp(2), "t_max": np.exp(3), "color": "#f39c12"},
    {"label": "ln(t) 3–4", "t_min": np.exp(3), "t_max": np.exp(4), "color": "#c0392b"},
    {"label": "ln(t) 4–5", "t_min": np.exp(4), "t_max": np.exp(5), "color": "#e74c3c"},
    {"label": "ln(t) 5–6", "t_min": np.exp(5), "t_max": np.exp(6), "color": "#2980b9"},
    {"label": "ln(t) 6–7", "t_min": np.exp(6), "t_max": np.exp(7), "color": "#27ae60"},
]

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

# --- 2. Plotting ---
results_df = pd.DataFrame(all_event_metrics)

if results_df.empty:
    print("No events detected. Check threshold parameters.")
else:
    wait_bins = np.geomspace(WAIT_MIN, WAIT_MAX, N_BINS + 1)

    # --- Figure 1: Individual subplots ---
    n_windows = len(WINDOWS)
    fig = make_subplots(
        rows=1, cols=n_windows,
        subplot_titles=[w["label"] for w in WINDOWS]
    )

    def power_law(x, scale, exponent):
        return scale * x ** exponent

    for col_idx, win in enumerate(WINDOWS, start=1):
        mask = (results_df["time"] >= win["t_min"]) & (results_df["time"] < win["t_max"])
        win_data = results_df[mask]["waiting_time"].dropna()

        if len(win_data) < 5:
            print(f"Not enough data in window {win['label']}")
            continue

        counts, _ = np.histogram(win_data, bins=wait_bins, density=True)

        fig.add_trace(
            go.Bar(
                x=wait_bins[:-1],
                y=counts,
                width=np.diff(wait_bins),
                marker_color=win["color"],
                opacity=0.7,
                name=win["label"],
                showlegend=False,
                hovertemplate="Wait: %{x:.4f}s<br>PDF: %{y:.4e}<extra></extra>"
            ),
            row=1, col=col_idx
        )

        x_fit = wait_bins[:-1][counts > 0]
        y_fit = counts[counts > 0]

        try:
            coeffs = np.polyfit(np.log(x_fit), np.log(y_fit), 1)
            fitted_exp   = coeffs[0]
            fitted_scale = np.exp(coeffs[1])

            print(f"{win['label']}: exponent = {fitted_exp:.3f}")

            x_line = np.geomspace(WAIT_MIN, WAIT_MAX, 200)
            y_line = power_law(x_line, fitted_scale, fitted_exp)

            fig.add_trace(
                go.Scatter(
                    x=x_line,
                    y=y_line,
                    mode='lines',
                    line=dict(color='black', width=2, dash='dash'),
                    showlegend=False,
                    hovertemplate=f"fit (α={fitted_exp:.2f}): %{{y:.4e}}<extra></extra>"
                ),
                row=1, col=col_idx
            )

            xref = "x domain" if col_idx == 1 else f"x{col_idx} domain"
            yref = "y domain" if col_idx == 1 else f"y{col_idx} domain"

            fig.add_annotation(
                x=0.05, y=0.05,
                xref=xref, yref=yref,
                text=f"α = {fitted_exp:.2f}",
                showarrow=False,
                font=dict(size=14, color="black"),
                bgcolor="white",
                bordercolor="black",
            )

        except Exception as e:
            print(f"Fit failed for window {win['label']}: {e}")

        fig.update_xaxes(title_text="Waiting Time (s)", type="log", row=1, col=col_idx)
        fig.update_yaxes(title_text="PDF (1/s)", type="log", row=1, col=col_idx)

    fig.update_layout(
        title_text="Waiting Time PDF by Time Window — Power Law Fits",
        showlegend=False,
        height=500,
        width=2300,
        template="plotly_white"
    )

    fig.show()

    # --- Figure 2: All windows overlaid with smoothed curves + fit lines ---
    SMOOTH_WINDOW = 5  # rolling average window in log space — increase for smoother

    fig2 = go.Figure()

    for win in WINDOWS:
        mask = (results_df["time"] >= win["t_min"]) & (results_df["time"] < win["t_max"])
        win_data = results_df[mask]["waiting_time"].dropna()
        if len(win_data) < 5:
            continue

        counts, _ = np.histogram(win_data, bins=wait_bins, density=True)
        x_fit = wait_bins[:-1][counts > 0]
        y_fit = counts[counts > 0]

        try:
            # Smooth in log space
            log_y_smooth = pd.Series(np.log(y_fit)).rolling(
                SMOOTH_WINDOW, center=True, min_periods=1
            ).mean()
            y_smooth = np.exp(log_y_smooth)

            # Smoothed data curve
            fig2.add_trace(go.Scatter(
                x=x_fit,
                y=y_smooth,
                mode='lines',
                line=dict(color=win["color"], width=2),
                name=win["label"],
                legendgroup=win["label"],
                showlegend=True,
            ))

            # Power law fit line (dashed, same color)
            coeffs = np.polyfit(np.log(x_fit), np.log(y_fit), 1)
            fitted_exp   = coeffs[0]
            fitted_scale = np.exp(coeffs[1])

            x_line = np.geomspace(x_fit.min(), x_fit.max(), 200)
            y_line = fitted_scale * x_line ** fitted_exp

            fig2.add_trace(go.Scatter(
                x=x_line,
                y=y_line,
                mode='lines',
                line=dict(color=win["color"], width=2, dash='dash'),
                name=f"{win['label']} α={fitted_exp:.2f}",
                legendgroup=win["label"],
                showlegend=True,
            ))

        except Exception as e:
            print(f"Overlay fit failed for {win['label']}: {e}")

    fig2.update_xaxes(title_text="Waiting Time (s)", type="log")
    fig2.update_yaxes(title_text="PDF (1/s)", type="log")
    fig2.update_layout(
        title_text="Power Law Fits — All Windows Overlaid",
        height=600,
        width=900,
        template="plotly_white",
        legend=dict(
            bordercolor="black",
            borderwidth=1,
        )
    )

    fig2.show()
    # fig.write_html("waiting_time_pdf_windows.html")
    # fig2.write_html("waiting_time_overlaid.html")