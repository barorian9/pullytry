import pandas as pd
import numpy as np
import plotly.graph_objects as go

rows = []
for i in range(1, 31):
    df = pd.read_csv(f"{i}.csv", names=["t_s","count","theta_rad","x_m"],
                     comment="#", on_bad_lines="skip", dtype=str)
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna()
    if df.empty:
        continue
    rows.append({
        "file_id": i,
        "x_total": abs(df["x_m"].iloc[-1] - df["x_m"].iloc[0]),
        "t_total": df["t_s"].iloc[-1] - df["t_s"].iloc[0],
    })

s = pd.DataFrame(rows)
print(s.describe())

fig = go.Figure(go.Scatter(x=s["t_total"], y=s["x_total"], mode="markers+text",
                           text=s["file_id"], textposition="top center"))
fig.update_layout(xaxis_title="run duration [s]", yaxis_title="total displacement [m]",
                  xaxis_type="log")
fig.show()

import pandas as pd
import numpy as np

df = pd.read_csv("event_data_shifted.csv")

# t=0 הוא האירוע הראשון בכל ריצה — לוג לא מוגדר שם
df = df[df["time"] > 0]
df["k"] = np.floor(np.log(df["time"])).astype(int)

summary = df.groupby("k").agg(
    n_events=("time", "count"),
    n_waits=("waiting_time", "count"),
    n_runs=("file_id", "nunique"),
)
summary["t_min"] = np.exp(summary.index)
summary["t_max"] = np.exp(summary.index + 1)
print(summary)

print("\nמשך כל ריצה (זמן מוסט של האירוע האחרון):")
print(df.groupby("file_id")["time"].max().describe())