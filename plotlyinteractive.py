import pandas as pd
import plotly.graph_objects as go

CSV = "7.2HZ.csv"

# Read CSV (skip junk lines, keep only valid numeric rows)
df = pd.read_csv(
    CSV,
    names=["t_s","count","theta_rad","x_m"],   # because your logger writes no header
    comment="#",
    on_bad_lines="skip"
)

# Convert to numeric and drop non-numeric rows
for c in df.columns:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna()

fig = go.Figure()
fig.add_trace(go.Scatter(
    x=df["t_s"],
    y=df["x_m"],
    mode="markers",      # use "markers" if you want every sample as dots
    name="x(t)"
))

fig.update_layout(
    title="Pulley displacement (raw)",
    xaxis_title="time [s]",
    yaxis_title="displacement x [m]",
    hovermode="x unified"
)




fig.show()
fig.write_html("pulley_displacement_interactive.html")
