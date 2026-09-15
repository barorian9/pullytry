> **SUPERSEDED.** This file predates `decisions_ledger_v2.md`. Claims in it
> about results, thresholds, aging, or Omori behaviour are not current. The
> ledger is authoritative. This file is kept for its procedural and code
> content only.

# Pulley Stick-Slip Analysis

Experimental data analysis testing whether a pulley under constant load, undergoing
stick-slip dynamics, shows self-organized criticality (SOC) and scale-free avalanche
statistics. Modeled after Shohat, Friedman & Lahini (2023) on crumpled sheets.

**Central question:** does a stick-slip pulley belong to the same universality class
as crumpled sheets?

---

## Documentation

| File | Contents |
|---|---|
| `project_summary.md` | **Master reference.** Methodology, parameters, every graph, results, physical interpretation, comparison to Shohat, current status. |
| `pipeline.txt` | Stage-by-stage description of the processing chain, plus the file dependency graph. |
| `README.md` | This index. |

---

## Data

| File | Contents |
|---|---|
| `1.csv` … `30.csv` | Raw runs. Columns: `t_s`, `count`, `theta_rad`, `x_m`. Lines starting with `#` are sync/comment lines. |
| `event_data.csv` | Detected slip events. Columns: `file_id`, `time`, `waiting_time`, `magnitude`. ~1572 events across 30 runs. **Working dataset for downstream analysis.** |
| `event_data_shifted.csv` | Same, with $t=0$ set to the first detected event of each run. Used for aging plots. |

Other CSVs (`3.5HZ.csv`, `7.2HZ.csv`, `noisecheck.csv`, `pulley_100Hz_long*.csv`,
`try1.csv`, `31.csv`) are single-run test captures from system characterization,
not part of the 30-run dataset.

---

## Code

### Pipeline

| Script | Reads | Writes | Role |
|---|---|---|---|
| `pulley_activate.py` | serial port | `N.csv` | Acquisition from Arduino encoder |
| `plotmultibple.py` | `1..30.csv` | `event_data.csv`, `event_data_shifted.csv` | **Event detection over all runs** + 3-panel plot (aging, $P(\Delta t)$, $P(M_x)$) |
| `log_for_dif_times.py` | `1..30.csv` | — | Aging analysis across five $\ln t$ windows. Re-runs detection independently. |
| `pm_finder.py` | `event_data.csv` | — | $P(M)$ with fixed 10 s threshold, plus count-per-size bar chart. **Current result.** |
| `inter_avalanche_claude.py` | `event_data_shifted.csv` | — | $P(M)$ with dynamic $C \cdot t$ threshold. Exploratory. |

### Diagnostics

| Script | Purpose |
|---|---|
| `cutoff_finder.py` | Per-run $\Delta t$ vs $t$ grid (30 panels) with fixed threshold candidates overlaid |
| `verifyeventdetection.py` | Single run: detected events overlaid on displacement and filtered velocity |
| `filtration.py` | Single run: raw vs median-filtered displacement |
| `filterossilations.py` | Single run: filter + detector, with merged events highlighted |

### Exploratory

Early characterization scripts, written before the pipeline was fixed. They use
different thresholds and most apply no ringing filter — **not comparable to pipeline
output**, kept for reference:

`derivative_threshold.py`, `timetry.py`, `tryploting.py`, `plotlyinteractive.py`, `12.py`

---

## Parameters

```python
VELOCITY_THRESHOLD  = 0.0085   # m/s  — minimum |v| to count a sample as active
COALESCENCE_TIME    = 0.135    # s    — merge active samples closer than this
RINGING_FREQ_HZ     = 6        # Hz   — median filter window = 1 period
FILE_COUNT          = 30

WAIT_MIN, WAIT_MAX  = 1e-3, 100    # s, log bin range
MAG_MIN,  MAG_MAX   = 1e-6, 1e-1   # m, log bin range
N_BINS              = 50

AVALANCHE_THRESHOLD = 10.0     # s    — fixed, not dynamic
```

These are duplicated across scripts rather than shared from one module. If you change
one, change all of them.

---

## Conventions

- Log-spaced bins (`np.geomspace`) for continuous distributions; integer bins for discrete $M$
- Histograms normalized to PDF (`density=True`)
- Power law fits via `np.polyfit` on $\log x$, $\log y$ — **never `curve_fit` in linear space**
- Plots in Plotly (`graph_objects`)

---

## Results so far

| Quantity | Result | Shohat (crumpled sheets) |
|---|---|---|
| Waiting time $P(\Delta t)$ | $\sim \Delta t^{-1}$ (Omori) | power law |
| Aging in $P(\Delta t)$ | none — $\alpha \approx -1.1$ stationary | present, $\alpha$ shifts |
| $\Delta t$ cutoff growth | $\sim t^{0.98}$ | $\sim t$ |
| Avalanche $P(M)$ | $\sim M^{-2.07}$ | $-2.1$ exp / $-1.8$ sim |

**Open issue:** sparse statistics on large avalanches — only 293 avalanches with
$M \geq 2$, and very few above $M = 10$. More runs needed before $\alpha$ is reliable.

---

## Reference

Shohat, Friedman & Lahini (2023), *Logarithmic aging via instability cascades in
disordered systems*, arXiv:2306.00567.
