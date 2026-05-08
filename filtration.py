import pandas as pd
import plotly.graph_objects as go
from scipy.signal import medfilt

# --- Parameters ---
CSV = "24.csv"
RINGING_FREQ_HZ = 6 # Change to 7.2 if that's your knocking frequency

# --- Load ---
df = pd.read_csv(CSV, names=["t_s", "count", "theta_rad", "x_m"], comment="#", on_bad_lines="skip")
for c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna()

# --- Filtration ---
dt = df["t_s"].diff().median()
fs = 1 / dt
window = int((1/RINGING_FREQ_HZ) * fs)
if window % 2 == 0: window += 1

# Apply Filter
df["x_filtered"] = medfilt(df["x_m"], kernel_size=window)

# --- Plotting ---
fig = go.Figure()

# 1. RAW Displacement (The "Before")
# Using a light color and thin line to show the ringing "cloud"
fig.add_trace(go.Scatter(
    x=df["t_s"], y=df["x_m"],
    name="Raw Displacement (with Ringing)",
    line=dict(color='red', width=1.5)
))

# 2. FILTERED Displacement (The "After")
# A bold dark line to show the physical path of the grains
fig.add_trace(go.Scatter(
    x=df["t_s"], y=df["x_filtered"],
    name="Filtered Displacement (Physics)",
    line=dict(color='black', width=2.5)
))

fig.update_layout(
    title=f"Displacement Filtration Check ({RINGING_FREQ_HZ}Hz Filter)",
    xaxis_title="Time [s]",
    yaxis_title="Position [m]",
    template="plotly_white",
    hovermode="x"
)

fig.show()