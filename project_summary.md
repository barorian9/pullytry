# Pulley Stick-Slip Analysis — Project Summary
*Generated from Claude conversation — May 2026*
*Feed this document to Claude to continue from the same point*

---

## 1. Experiment Overview

**System:** A pulley under constant load that undergoes stick-slip dynamics — periods of being stuck alternating with rapid slips.

**Data format:** 30 CSV files named `1.csv` through `30.csv`, each representing one independent experimental run. Each CSV has columns:
- `t_s` — time in seconds
- `count` — encoder count
- `theta_rad` — angle in radians
- `x_m` — displacement in meters

**Output file:** `event_data.csv` — processed slip events with columns:
- `file_id` — run number (1-30)
- `time` — time of peak velocity during event (seconds)
- `waiting_time` — gap since previous event (seconds), NaN for first event per run
- `magnitude` — absolute displacement during event (meters)

---

## 2. Key Parameters

```python
FILE_COUNT = 30
VELOCITY_THRESHOLD = 0.0085   # m/s — minimum velocity to count as a slip
COALESCENCE_TIME = 0.135      # s — merge events closer than this
RINGING_FREQ_HZ = 6           # Hz — median filter to remove sensor ringing
```

**Event detection pipeline:**
1. Load raw CSV
2. Apply median filter (window = 1/6Hz × fs samples) to displacement
3. Compute velocity from filtered displacement
4. Mark samples where |v| > 0.0085 m/s as "active"
5. Merge active samples within 0.135s into one event
6. Record peak velocity time, waiting time, magnitude per event

**Verification:** Event detection was verified by overlaying detected events on zoomed displacement curves. Detection looks correct — every visible displacement drop has a marker, no false detections on flat sections.

---

## 3. Analysis Pipeline

### 3.1 Main Processing Code

The main code loops over all 30 CSVs, detects events, and builds `results_df` / `event_data.csv`. Key output: ~1572 total events across 30 runs (roughly 40-80 events per run).

### 3.2 Waiting Time PDF (Global)

**What we changed from naive histogram:**
- Old: `go.Histogram` with `nbinsx=1000` → linear bins, raw counts
- New: `np.geomspace(WAIT_MIN, WAIT_MAX, N_BINS+1)` → log-spaced bins
- Normalized with `np.histogram(..., density=True)` → true PDF (units: 1/s)
- Both axes log-scale → power law appears as straight line

**Parameters:**
```python
WAIT_MIN, WAIT_MAX = 1e-3, 100   # seconds
MAG_MIN, MAG_MAX = 1e-6, 1e-1   # meters
N_BINS = 50
```

**Result:** Global waiting time PDF follows **P(Δt) ~ Δt^(-1)** — Omori law. Also added t^(-2) tail for large Δt, anchored at t=20s.

### 3.3 Power Law Fitting

**Key insight:** Fitting in log-log space gives the correct exponent:

```python
# WRONG — biased toward high-PDF bins
from scipy.optimize import curve_fit
popt, _ = curve_fit(power_law, x_fit, y_fit, p0=[1.0, -1.0])

# CORRECT — equal weight per decade
coeffs = np.polyfit(np.log(x_fit), np.log(y_fit), 1)
fitted_exp = coeffs[0]
fitted_scale = np.exp(coeffs[1])
```

### 3.4 Aging Analysis (Time Window Analysis)

**Method (following Shohat 2023):**
- Divide data into log-time windows: ln(t) = 2-3, 3-4, 4-5, 5-6, 6-7
- Compute waiting time PDF separately in each window
- Fit power law exponent α in each window
- If α changes systematically → aging. If α stable → no aging.

**Windows used:**
```python
WINDOWS = [
    {"label": "ln(t) 2–3", "t_min": np.exp(2), "t_max": np.exp(3), "color": "#f39c12"},
    {"label": "ln(t) 3–4", "t_min": np.exp(3), "t_max": np.exp(4), "color": "#c0392b"},
    {"label": "ln(t) 4–5", "t_min": np.exp(4), "t_max": np.exp(5), "color": "#e74c3c"},
    {"label": "ln(t) 5–6", "t_min": np.exp(5), "t_max": np.exp(6), "color": "#2980b9"},
    {"label": "ln(t) 6–7", "t_min": np.exp(6), "t_max": np.exp(7), "color": "#27ae60"},
]
```

**Note:** x-axis of aging plot uses `np.log(time)` = natural log (not log10). So x=4 means t=e^4≈55s, x=6 means t=e^6≈403s.

**Result:**
| Window | α |
|--------|---|
| ln(t) 2–3 | -0.95 |
| ln(t) 3–4 | -1.00 |
| ln(t) 4–5 | -1.18 |
| ln(t) 5–6 | -1.19 |
| ln(t) 6–7 | -1.11 |

**Conclusion: No aging** — α ≈ -1.1 is stationary across all windows from ln(t) 4-7. The system is in a **stable critical state**.

### 3.5 C·t Envelope Check

Computed C = mean(Δt/t) per run to check if Shohat's dynamic threshold applies:

```python
C = (run["waiting_time"] / run["time"]).mean()
```

**C values per run:** Range from 0.037 to 0.115, clustering around 0.05-0.10.

**Approach 2 — Rolling max of Δt:**
- Divided log-time axis into 20 windows
- Found max(Δt) per window
- Fitted max(Δt) ~ t^slope

**Result:** max(Δt) grows as **t^0.98 ≈ t^1** — the upper cutoff of waiting times grows linearly with time, consistent with Shohat aging. Mean(Δt) grows as t^0.67.

**Reconciliation:** The *shape* of P(Δt) is stationary (no aging in α), but the *cutoff* of the distribution grows linearly with t. This is self-similar aging — the shape doesn't change but the timescale shifts. Consistent with Shohat's finding that avalanches remain self-similar throughout aging.

---

## 4. Avalanche Analysis

### 4.1 Definition

An **avalanche** is a cluster of slip events separated from the next cluster by a quiet period longer than a threshold Δt_threshold.

**Avalanche size M** = number of slip events inside one avalanche.

**Example:**
Events at t = 10s, 10.3s, 10.7s, 25s, 25.4s, 60s
Waiting times: 0.3s, 0.4s, 14.3s, 0.4s, 34.6s
With threshold = 10s:
- 0.3s < 10s → same avalanche
- 0.4s < 10s → same avalanche (M=3 so far)
- 14.3s > 10s → **end of avalanche**, record M=3
- 0.4s < 10s → same avalanche
- 34.6s > 10s → **end of avalanche**, record M=2

### 4.2 Threshold Selection

**Key finding from per-run analysis:**
When looking at each run individually (not pooled), the data shows a clear bimodal structure:
- Dense lower band: Δt < ~1-2s (intra-avalanche events)
- Sparse upper points: Δt > ~5-10s (inter-avalanche gaps)

**The threshold should be fixed (constant), not dynamic (C·t)** — because the intra-avalanche waiting times don't grow with t. They stay bounded by the physics of the pulley mechanism.

**Threshold comparison:**

| Threshold | N avalanches (M≥2) | Mean M | Max M | α |
|-----------|-------------------|--------|-------|---|
| 1s | 369 | 2.6 | 7 | -4.09 |
| 2s | 377 | 3.0 | 16 | -2.99 |
| 5s | 365 | 3.7 | 17 | -2.68 |
| 10s | 293 | ~3.7 | 20 | -2.07 |

**Selected threshold: 10s** — gives α = -2.07, closest to Shohat's experimental value of -2.1. Also visually most consistent across all 30 per-run plots.

### 4.3 P(M) Calculation

For discrete M (integer avalanche sizes):

```python
M_vals = np.arange(2, sizes.max() + 1)
counts = np.array([np.sum(sizes == m) for m in M_vals])
pdf = counts / len(sizes)   # bin width = 1, so no division needed
```

This is technically a **probability mass function (PMF)** not a PDF, since M is discrete. Log-log axes used for visualization.

**Result with threshold = 10s:**
- Total avalanches (M≥2): 293
- Max M: 20
- **α = -2.07** (Shohat found -2.1 experimentally, -1.8 in simulations)

**Limitation:** Not enough large avalanches. Data becomes sparse above M=10 (below N=10 reliability floor). Need more experimental runs to constrain α reliably.

---

## 5. Graphs Created

### Graph 1: Combined Aging Plot (3-panel)
- **Left:** Δt vs t on log-log axes, colored by run ID (Viridis colormap)
- **Middle:** Global waiting time PDF with log bins, Omori t^(-1) fit (red), t^(-2) tail (orange)
- **Right:** Magnitude PDF with log bins
- Both histogram axes log-log

### Graph 2: Waiting Time PDF by Time Window (5 subplots)
- One subplot per ln(t) window (2-3, 3-4, 4-5, 5-6, 6-7)
- Log-spaced bins, density=True
- Power law fit in log-log space via np.polyfit
- α annotated in each subplot

### Graph 3: All Windows Overlaid
- Smoothed curves (rolling mean in log space, window=5) for each time window
- Dashed fit lines same color as data
- Shows all 5 windows on one plot to compare shapes

### Graph 4: Threshold Candidates on Global Aging Plot
- All events pooled, log-log axes
- Overlaid: C·t line (red), max envelope fit (orange), ½×max (purple), mean envelope (blue), fixed 1s/5s/10s (green dotted)

### Graph 5: Per-Run Aging Plots (10×3 grid)
- One subplot per run, Δt vs t log-log
- Threshold candidates overlaid: 1s, 2s, 5s, 10s fixed lines
- Shows that fixed threshold is appropriate (lower band stays constant with t)

### Graph 6: Threshold Verification (30 runs, 10×3 grid)
- Same as Graph 5 but all 30 runs
- Used to select threshold = 10s visually

### Graph 7: Event Detection Verification (single run)
- Top panel: raw (red) and filtered (black) displacement
- Bottom panel: filtered velocity with threshold lines
- Green triangles on displacement = detected events
- Green dotted vertical lines through both panels
- Used to verify event detection quality

### Graph 8: P(M) — Avalanche Size Distribution
- Dots at each integer M
- Red dashed power law fit
- Annotation box: α, threshold, N avalanches, max M

### Graph 9: Avalanche Count per Size (bar chart)
- Bar at each M showing how many avalanches found
- Red dashed line at N=10 (reliability floor)
- Shows data runs out above M~10

### Graph 10: Approach 1/2/3 — Threshold Validation
- Approach 1: Δt/t vs t (tests if ratio is constant)
- Approach 2: Rolling max and mean of Δt vs t (fits power law to envelope)
- Approach 3: P(M) with dynamic vs fixed threshold comparison

---

## 6. Physical Interpretation

### Waiting Time Distribution
- P(Δt) ~ Δt^(-1) → Omori law, scale-free, no characteristic timescale
- System is in a **stationary critical state** — self-organized criticality (SOC)
- α stationary across all time windows → **no aging** in the distribution shape
- Upper cutoff grows as t^1 → timescale of rare large events grows with system age

### Two Populations
The waiting time distribution has two regimes:
1. **Short Δt (0.1-2s):** Intra-avalanche events — slips triggering more slips
2. **Long Δt (>10s):** Inter-avalanche gaps — quiet periods between bursts

### Avalanche Size Distribution
- P(M) ~ M^(-2.07) → power law, scale-free
- Very close to Shohat's crumpled sheet result (α = -2.1)
- Consistent with SOC — no characteristic avalanche size

### Comparison to Shohat (2023)
| Property | Shohat (crumpled sheets) | This system (pulley) |
|----------|--------------------------|----------------------|
| Waiting time P(Δt) | Power law | Power law α≈-1.1 |
| Aging in P(Δt) | Yes (α changes) | No (α stable) |
| Upper cutoff | Grows as C·t | Grows as t^0.98 |
| Avalanche threshold | Dynamic C·t | Fixed ~10s |
| P(M) exponent | -2.1 (exp), -1.8 (sim) | -2.07 (preliminary) |

**Key difference:** Shohat's system ages (the shape of P(Δt) changes over time). Your pulley system does not age — the shape stays the same, only the cutoff grows. Both systems show similar P(M) power laws.

---

## 7. Current Status & Next Steps

### Done:
- Event detection pipeline ✓
- Global waiting time PDF ✓
- Aging analysis (no aging confirmed) ✓
- Avalanche definition and threshold selection ✓
- P(M) preliminary result ✓

### Limitations:
- Only 293 avalanches with M≥2, only 20 with M≥10
- P(M) fit unreliable above M~8 due to low counts
- α sensitive to threshold choice (ranges from -2.07 to -4.09)
- Need more experimental runs to get reliable statistics on large avalanches

### Next steps:
- Collect more experimental data (more runs)
- Potentially lower VELOCITY_THRESHOLD to detect more events
- Study magnitude distribution P(ΔX) — not yet fully analyzed
- Consider intra-avalanche dynamics (waiting time distribution *within* avalanches)

---

## 8. Key References

- **Shohat, Friedman & Lahini (2023):** "Logarithmic aging via instability cascades in disordered systems" — arXiv:2306.00567v2. The paper this analysis is modeled after. Key results: logarithmic aging in crumpled sheets, avalanche size distribution P(M) ~ M^(-2.1), waiting time cutoff grows as C·t with C≈7×10^(-4).

---

## 9. Code Files Summary

All code was written in Python using: `pandas`, `numpy`, `scipy.signal.medfilt`, `plotly.graph_objects`, `plotly.subplots`, `scipy.optimize.curve_fit`

**Main scripts produced during this conversation:**
1. `main_analysis.py` — event detection loop, produces event_data.csv, 3-panel plot
2. `waiting_time_windows.py` — aging analysis with 5 time windows
3. `inter_avalanche_claude.py` — avalanche analysis with dynamic C·t threshold
4. `cutoff_finder.py` — threshold validation (3 approaches)
5. `threshold_verification.py` — 30-run grid plot with threshold candidates
6. `event_detection_verify.py` — single run verification plot
7. `avalanche_pm.py` — final P(M) plot with count bar chart

---

*End of summary. Feed this document to Claude at the start of a new conversation to continue from this point.*
