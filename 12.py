import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("pulley_100Hz_long.csv")

plt.figure(figsize=(10,4))
plt.plot(df["t_s"], df["x_m"], linewidth=0.8)
plt.xlabel("Time [s]")
plt.ylabel("Displacement x [m]")
plt.grid()
plt.tight_layout()
()
plt.show()
t0 = 10     # start time [s]
t1 = 12     # end time [s]

zoom = df[(df["t_s"] > t0) & (df["t_s"] < t1)]

plt.figure(figsize=(10,4))
plt.plot(zoom["t_s"], zoom["x_m"], marker='.', linestyle='-')
plt.xlabel("Time [s]")
plt.ylabel("Displacement x [m]")
plt.grid()
plt.tight_layout()
plt.show()
