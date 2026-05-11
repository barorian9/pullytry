import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

df = pd.read_csv("event_data.csv")
df = df.dropna(subset=["waiting_time", "time"])
df = df[df["waiting_time"] > 0]
df = df[df["time"] > 0]

THRESHOLD = 10.0

# --- Get avalanche sizes ---
all_sizes = []
for run_id, run in df.groupby("file_id"):
    run = run.sort_values("time").reset_index(drop=True)
    current = 1
    for i in range(1, len(run)):
        if run.loc[i, "waiting_time"] > THRESHOLD:
            all_sizes.append(current)
            current = 1
        else:
            current += 1
    all_sizes.append(current)

all_sizes = np.array(all_sizes)
sizes = all_sizes[all_sizes >= 2]

print(f"Total avalanches (M>=2): {len(sizes)}")
print(f"Max M: {sizes.max()}")
for m in range(2, all_sizes.max() + 1):
    c = np.sum(all_sizes == m)
    if c > 0:
        print(f"  M={m}: {c}")

# PDF
M_vals = np.arange(2, sizes.max() + 1)
counts = np.array([np.sum(sizes == m) for m in M_vals])
pdf = counts / len(sizes)
mask = counts > 0
x_fit = M_vals[mask].astype(float)
y_fit = pdf[mask]

# Power law fit
coeffs = np.polyfit(np.log(x_fit), np.log(y_fit), 1)
alpha = coeffs[0]
scale = np.exp(coeffs[1])
x_line = np.geomspace(2, sizes.max(), 200)
y_line = scale * x_line ** alpha

# ============================================================
# Figure 1: Count per avalanche size — shows data limitation
# ============================================================
fig1 = go.Figure()

fig1.add_trace(go.Bar(
    x=M_vals[mask],
    y=counts[mask],
    marker_color="#2980b9",
    text=counts[mask],
    textposition="outside",
    showlegend=False
))

fig1.add_hline(
    y=10,
    line=dict(color="red", width=1.5, dash="dash"),
    annotation_text="N = 10",
    annotation_position="top right"
)

fig1.update_xaxes(
    title_text="Avalanche size M (events)",
    tickvals=list(M_vals[mask])
)
fig1.update_yaxes(title_text="Count")
fig1.update_layout(
    title_text=f"Avalanche count per size — threshold {THRESHOLD}s | 30 runs | {len(sizes)} avalanches total",
    height=450, width=800,
    template="plotly_white"
)
fig1.show()

# ============================================================
# Figure 2: P(M) — simple and clean
# ============================================================
fig2 = go.Figure()

fig2.add_trace(go.Scatter(
    x=x_fit,
    y=y_fit,
    mode="markers",
    marker=dict(color="#2c3e50", size=9),
    name="P(M)",
    hovertemplate="M=%{x}<br>P(M)=%{y:.4f}<extra></extra>"
))

fig2.add_trace(go.Scatter(
    x=x_line,
    y=y_line,
    mode="lines",
    line=dict(color="#e74c3c", width=2, dash="dash"),
    name=f"fit: α = {alpha:.2f}"
))

fig2.add_annotation(
    x=0.97, y=0.97,
    xref="paper", yref="paper",
    text=f"α = {alpha:.2f}<br>threshold = {THRESHOLD}s<br>N avalanches = {len(sizes)}<br>Max M = {sizes.max()}",
    showarrow=False,
    align="right",
    bgcolor="white",
    bordercolor="black",
    borderwidth=1,
    font=dict(size=12)
)

fig2.update_xaxes(
    title_text="Avalanche size M (events)",
    type="log",
    tickvals=[2, 3, 4, 5, 6, 7, 8, 10, 12, 15]
)
fig2.update_yaxes(title_text="PDF P(M)", type="log")
fig2.update_layout(
    title_text="Avalanche size distribution P(M) — preliminary",
    height=500, width=750,
    template="plotly_white"
)

fig2.update_yaxes(
    title_text="PDF P(M)",
    type="log",
    range=[-3, 0]  # מ-0.001 עד 1
)
fig2.show()