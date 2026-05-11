import pandas as pd
import numpy as np
import plotly.graph_objects as go

df = pd.read_csv("event_data_shifted.csv")
df = df.dropna(subset=["waiting_time", "time"])
df = df[df["waiting_time"] > 0]
df = df[df["time"] > 0]

all_avalanche_sizes = []

for run_id, run in df.groupby("file_id"):
    run = run.sort_values("time").reset_index(drop=True)

    # Compute C for this run
    C = (run["waiting_time"] / run["time"]).mean()
    print(f"Run {run_id}: C = {C:.4f}, N events = {len(run)}")

    # Define avalanche threshold: Δt > C·t means new avalanche
    current_avalanche = 1
    for i in range(1, len(run)):
        dt = run.loc[i, "waiting_time"]
        t  = run.loc[i, "time"]
        threshold = C * t
        if dt > threshold:
            all_avalanche_sizes.append(current_avalanche)
            current_avalanche = 1
        else:
            current_avalanche += 1
    all_avalanche_sizes.append(current_avalanche)

sizes = np.array(all_avalanche_sizes)
print(f"\nTotal avalanches: {len(sizes)}")
print(f"Mean size: {sizes.mean():.2f}")
print(f"Max size: {sizes.max()}")

# Filter M >= 2
sizes = sizes[sizes >= 2]
print(f"Avalanches with M >= 2: {len(sizes)}")

# Integer bins (M is discrete)
M_vals = np.arange(2, sizes.max() + 1)
counts_raw = np.array([np.sum(sizes == m) for m in M_vals])
total = len(sizes)
pdf = counts_raw / total  # bin width = 1

mask = counts_raw > 0
x_fit = M_vals[mask].astype(float)
y_fit = pdf[mask]

# Power law fit in log-log space
coeffs = np.polyfit(np.log(x_fit), np.log(y_fit), 1)
alpha = coeffs[0]
scale = np.exp(coeffs[1])
print(f"\nPower law exponent: α = {alpha:.3f}")

x_line = np.geomspace(2, sizes.max(), 200)
y_line = scale * x_line ** alpha

# --- Plot ---
fig = go.Figure()

fig.add_trace(go.Scatter(
    x=x_fit,
    y=y_fit,
    mode='markers',
    marker=dict(color='#e74c3c', size=8),
    name="P(M)",
    hovertemplate="M: %{x}<br>PDF: %{y:.4e}<extra></extra>"
))

fig.add_trace(go.Scatter(
    x=x_line,
    y=y_line,
    mode='lines',
    line=dict(color='black', width=2, dash='dash'),
    name=f"fit: M^{alpha:.2f}",
))

fig.add_annotation(
    x=0.95, y=0.95,
    xref="paper", yref="paper",
    text=f"α = {alpha:.2f}",
    showarrow=False,
    font=dict(size=14),
    bgcolor="white",
    bordercolor="black",
    borderwidth=1,
    align="right"
)

fig.update_xaxes(
    title_text="Avalanche size M (events)",
    type="log",
    tickvals=[2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 40]
)
fig.update_yaxes(
    title_text="PDF P(M)",
    type="log"
)
fig.update_layout(
    title_text="Avalanche size distribution P(M) —  dynamic threshold",
    height=500,
    width=800,
    template="plotly_white",
    showlegend=True
)

fig.show()