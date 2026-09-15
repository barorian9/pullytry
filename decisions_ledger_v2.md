# Pulley Stick-Slip — Decisions Ledger v2

*Corrected 2026-09-14 by going through the v1 ledger item by item with Bar. This supersedes v1 and supersedes the "Key results" section of `project_summary.md` wherever they disagree.*

**Status legend:** `LOCKED` = decided and justified · `WORKING` = in use but not decided · `OPEN` = known unresolved · `NOT FOUND` = never established, do not reconstruct

---

## A. Locked decisions

### A1. PDF normalization — `LOCKED`

| | |
|---|---|
| **Continuous** (`waiting_time`, `magnitude`) | `np.histogram(..., density=True)` → count / (total × bin's own linear width). Units 1/s and 1/m. |
| **Discrete** (avalanche size $M$) | `counts / len(sizes)`, bin width 1 → a PMF, not a PDF. |
| **Why** | $y$-axis independent of bin count and sample size; $\int P\,dx = 1$. Log bins have unequal *linear* widths, so per-bin division is required — `density=True` does this automatically. For integer $M$ the width is 1, so the division is a no-op. |

### A2. Log binning — `LOCKED`

- `np.geomspace(MIN, MAX, N_BINS+1)`, `N_BINS = 50`
- Waiting time: `WAIT_MIN, WAIT_MAX = 1e-3, 100` s
- Magnitude: `MAG_MIN, MAG_MAX = 1e-6, 1e-1` m
- Empty bins excluded from fits via `counts > 0`; left as gaps in the plot ($\log 0$ undefined)
- **Discrete $M$: integer bins** `np.arange(2, sizes.max()+1)`, never geomspace

**Revision on record:** log-binning $M$ gave $\alpha = -2.40$; integer binning gave $-1.61$ on identical data. Integer kept as the honest choice for discrete counts.

**⚠ Flagged against A6:** `WAIT_MIN = 1e-3` s is below the physical floor. See D2.

### A3. Global $t_0$ — `LOCKED`

The logger clock starts at drive activation plus a manual keypress — **not uniform across runs**. Therefore a new origin was defined:

$$t_0^{(r)} = \min_{j} t_j^{(r)} \quad \text{(first detected event in run } r\text{)}, \qquad t_{\text{shifted}} = t - t_0^{(r)}$$

Written to `event_data_shifted.csv`; raw `event_data.csv` preserved.

| Analysis | Time base |
|---|---|
| Aging windows ($\ln t$) | shifted |
| Dynamic threshold $C \cdot t$ | shifted |
| $P(\Delta t)$, $P(\text{magnitude})$ | raw |
| $P(M)$, fixed 10 s | raw |

**Not an inconsistency:** $\Delta t$ and magnitude are origin-independent, so the raw/shifted split has no effect on them. The origin only matters for the aging analysis and $C\cdot t$, and both already use shifted time.

**Known weakness:** $t_0$ is itself a random draw from $P(\Delta t)$, and it depends on `VELOCITY_THRESHOLD` — lowering the threshold moves $t_0$ earlier in every run. The earliest windows ($\ln t = 2$–$3$, $3$–$4$) are correspondingly less trustworthy. **New rig: use an external timestamp of load application.**

### A4. Magnitude — filtered, deliberate — `LOCKED`

Two quantities share the letter $M$ — keep them apart:

| Symbol | Definition |
|---|---|
| `magnitude` column | $\lvert x_{\text{filtered}}(\text{end}) - x_{\text{filtered}}(\text{start})\rvert$, metres |
| $M$ in $P(M)$ | number of detected events inside one avalanche (dimensionless count), $M \ge 2$ |

**Why filtered and not raw — resolved:** after every slip the string rings — sinusoidal oscillations clearly visible in raw `x_m`, at the frequency of the drive. A median filter of width $\lfloor f_s / 6 \rfloor = 17$ samples ($0.17$ s $= 1/6$ Hz) removes them. Without it the ringing would register as events and $P(\Delta t)$ would show a spurious peak near $1/6$ s. **The filtering is necessary, not incidental.**

**Side effect, for the record:** a median filter suppresses fast excursions, so it systematically shrinks the largest and fastest slips — the tail of the magnitude distribution. Not a reason to remove it; a reason not to over-read $P(\Delta X)$.

`magnitude` has only ever been *plotted* (panel 3). It was never fitted. All $P(M)$ work uses event counts.

### A5. $M = 1$ discarded — `LOCKED`
Single events are not avalanches. Only $M \ge 2$ enters $P(M)$.

**A5b.** The first detected event of each run is excluded from avalanche construction (`dropna(subset=["waiting_time"])` in both avalanche scripts). It has no defined preceding gap — the record is left-truncated at $t_0$, which is that event itself. Bound on the effect: $\le 30$ of 293 avalanches shift by one, $\lvert\delta\alpha_M\rvert \lesssim 0.1$. Deliberate, not a bug.

### A6. Sampling rate — `LOCKED` *(measured 2026-09-14)*

$$f_s = 100\ \text{Hz, uniform across all 30 runs}$$

- $\Delta t_{\text{sample}}$: median $0.01000$ s, 1st percentile $0.00900$, 99th $0.01100$ — steady within runs (one outlier in run 4)
- Sampling is **clock-driven, not encoder-driven**
- Median filter window = **17 samples in every run** → no cross-run non-uniformity in detection

### A7. Fitting — `LOCKED`

```python
coeffs = np.polyfit(np.log(x), np.log(y), 1)
alpha  = coeffs[0]
scale  = np.exp(coeffs[1])
```

- **Revision:** `scipy.optimize.curve_fit` in linear space was used first and abandoned — it is dominated by large $y$ at small $x$ and effectively ignores the tail; the fitted line floated visibly above the histogram. Log-log least squares weights each decade equally.
- Fit range = all non-empty bins across the full binning range. **No cutoff exclusion applied.** A proposal to drop the first few waiting-time bins was not adopted; the sparse large-$M$ tail was retained in the $P(M)$ fit.

### A8. Runs treated as independent — `LOCKED`
Avalanches are built per `file_id`. Pooling happens only afterwards, at the level of the resulting $M$ values.

---

## B. Working defaults — in use, NOT decided

### B1. Avalanche threshold — `WORKING` / `OPEN`

**Current working choice: fixed $\Delta t_{\text{thr}} = 10$ s.** A gap $> 10$ s ends the avalanche.

This is a working default, **not a conclusion.** The lean toward a fixed rather than dynamic threshold rests on per-run inspection (the lower $\Delta t$ band does not grow with $t$) — real evidence, not decisive. The value 10 itself was never derived from any criterion.

**The full arc, for reference:**

| Stage | What happened |
|---|---|
| (a) | Avalanche threshold initially conflated with `COALESCENCE_TIME = 0.135` s — corrected; that parameter is event detection only |
| (b) | Shohat's dynamic $C\cdot t$ tried; $C = \langle \Delta t / t\rangle$ per run, range $0.037$–$0.115$; gave $\alpha = -2.40$ (log bins) → $-1.61$ (integer bins) |
| (c) | Validity challenged three ways: $\Delta t / t$ vs $t$; rolling max/mean per log-window ($\max\Delta t \sim t^{0.98}$, $\mathrm{mean}\,\Delta t \sim t^{0.67}$); dynamic vs fixed-5 s $P(M)$ |
| (d) | Candidates overlaid on pooled aging plot — no clean gap |
| (e) | Pooling identified as the problem (runs independent, $t_0$ jitter) → moved to per-run plots |
| (f) | Per-run: lower $\Delta t$ band does not grow with $t$ → lean toward fixed |
| (g) | Fixed candidates: 1 s → $\alpha = -4.09$ (max $M = 7$); 2 s → $-2.99$ (16); 5 s → $-2.68$ (17); **10 s → $-2.07$ (20, 293 avalanches)** |

**⚠ Do not cite:** the phrase *"10 s was chosen on the basis of the graphs"* appears only in a Hebrew message drafted for the supervisor. It is Claude's wording, accepted in passing — **not a justification that was ever reached.**

**⚠ Consequence for reporting:** $\alpha_M = -2.07$ is a function of an open parameter and cannot stand as a result. The correct object is the curve $\alpha_M(\Delta t_{\text{thr}})$, not a number. Note the monotone slide $-4.09 \to -2.99 \to -2.68 \to -2.07$.

**Code status:** both implementations live. `pm_finder.py` = fixed 10 s on raw. `inter_avalanche_claude.py` = dynamic $C\cdot t$ on shifted. The second is a **reference implementation, not dead code.**

**Suspected bug in the dynamic branch:** $C = \langle \Delta t / t\rangle$ is computed on *shifted* time, where the second event of each run sits at very small $t_{\text{shifted}}$ → its ratio blows up and inflates $C$. May explain the wide spread in $C$.

### B2. Detection parameters — values `WORKING`, provenance partly `NOT FOUND`

```python
VELOCITY_THRESHOLD = 0.0085   # m/s, on |v| from filtered displacement
COALESCENCE_TIME   = 0.135    # s
RINGING_FREQ_HZ    = 6        # Hz → median filter kernel floor(fs/6), forced odd
```

| Parameter | Provenance |
|---|---|
| `RINGING_FREQ_HZ = 6` | **Established (A4)** — measured string ringing frequency, visible in raw traces |
| `VELOCITY_THRESHOLD` | **NOT FOUND** |
| `COALESCENCE_TIME` | **NOT FOUND** |

Both arrived pre-set in the original code paste. What exists instead: a post-hoc **visual verification** — events overlaid on filtered displacement for run 1 at four zoom levels — after which all three were left unchanged (*"a few misses but in general it looks good"*).

An older script used `COALESCENCE_TIME = 0.065` and a negative-only test (`v < -threshold`); confirmed it did **not** produce `event_data.csv`.

---

## C. Measured facts *(new, 2026-09-14)*

### C1. Event counts and run coverage

Windows are in shifted time, $k$ such that $t \in [e^k, e^{k+1}]$:

| $k$ | $t$ range [s] | events | runs contributing |
|---|---|---|---|
| −2 | 0.14–0.37 | 6 | 6 |
| −1 | 0.37–1.0 | 13 | 12 |
| 0 | 1.0–2.7 | 14 | 9 |
| 1 | 2.7–7.4 | 19 | 10 |
| 2 | 7.4–20 | 42 | 18 |
| 3 | 20–55 | 99 | 19 |
| 4 | 55–148 | 275 | 24 |
| 5 | 148–403 | **644** | **30** |
| 6 | 403–1097 | 413 | 23 |
| 7 | 1097–2981 | 47 | 3 |

`n_waits` equals `n_events` in every window.

**Reading:** every run survives past 148 s (min duration $153$ s), so the missing runs in $k = 2,3,4$ are runs with **no events in that window**, not runs that ended. Real censoring appears only at the top: median run duration is $576$ s, so roughly half the runs die *inside* $k=6$. **$k = 5$ is the only window that is both fully populated and uncensored.**

### C2. Run durations and total displacement

| | |
|---|---|
| Duration | min $153$ s, median $576$ s, mean $668$ s, max $2815$ s — a factor of $\approx 19$ spread |
| Total displacement $x_{\text{total}}$ | $\approx 0.61$ m, 24 of 30 runs within $0.55$–$0.66$ m, **no trend against duration** |
| Outliers | runs 2 ($0.78$), 8 ($0.80$), 23 ($0.83$) — $\approx 30\%$ above median |

All 30 runs cover the **same distance** in times spanning a factor of 19. What varies between runs is not how much happened, but how slowly.

---

## D. Open issues

### D1. Data-dependent censoring — `OPEN`
Every run stops when the weight reaches the floor, i.e. at a cumulative-sum crossing $\sum_j M_j \approx L$ with $L \approx 0.61 \pm 0.05$ m (small variation attributed to non-uniform starting height). Since run duration varies by $19\times$ at fixed distance, **the late windows are populated by the subset of slow runs.** Direction of bias: toward apparent aging. **Not quantified.** Must be revisited before any time-dependence claim is published. *New rig: a track long enough that most runs end on the timer, not the floor — this removes the problem entirely.*

### D2. Encoder resolution vs velocity threshold — `OPEN`
$\texttt{VELOCITY\_THRESHOLD} = 0.0085$ m/s at $f_s = 100$ Hz $\Rightarrow$ $8.5\times10^{-5}$ m between consecutive samples. **Unknown how many encoder counts that is.** If 1–2, the threshold sits on the quantization floor rather than on physical noise, and any slip smaller than one count is invisible by construction — a hard lower bound on $P(M)$ and $P(\Delta X)$.
*To get it:* counts per revolution × pulley diameter → metres per count. Or from the data: smallest non-zero difference in `x_m` within a run.

### D3. $\max\Delta t \sim t^{0.98}$ may be a binning artifact — `OPEN`
$\max\Delta t$ is measured *inside* a window $[e^k, e^{k+1}]$ whose linear width is $e^k(e-1) \propto t$. No gap can exceed the window containing it, so the upper bound on $\max\Delta t$ grows linearly with $t$ **by construction**. With a heavy-tailed distribution the maximum tracks that bound and returns an exponent $\approx 1$. The measured value is $0.98$.
*Decisive test (cheap):* shuffle the $\Delta t$ values across times, keeping the time stamps, and recompute $\max\Delta t$ per window. The shuffled data has no aging by construction. If $t^{\approx 1}$ survives, the finding is dead. Run the same test on $\mathrm{mean}\,\Delta t \sim t^{0.67}$.

### D4. No uncertainty on any exponent — `OPEN`
Not computed anywhere: not for $\alpha$ of $P(\Delta t)$, not for $\alpha_M$, not for the $0.98$ exponent. **"No aging" is therefore not a result** — it is a visual observation on five numbers ($-0.95, -1.00, -1.18, -1.19, -1.11$), whose ordering is monotone rather than scattered. When it is addressed: bootstrap over raw events, not `polyfit` errors on bin points (bins are neither independent nor equally weighted).

The two shallowest values ($-0.95$, $-1.00$) come from the two least-supported windows: $k=2$ has 42 events from 18 runs, $k=3$ has 99 from 19 (C1), and A3 notes $t_0$ jitter hits early windows hardest. The apparent monotone trend rests on exactly those two.

### D5. Effective time resolution — `OPEN` *(minor)*
The median filter is $0.17$ s wide and `COALESCENCE_TIME` is $0.135$ s, so two slips separated by less than $\approx 0.17$ s cannot be resolved. **The true lower limit of $P(\Delta t)$ is $\approx 0.17$ s, not `WAIT_MIN = 1e-3` s.** Bins below $0.01$ s are empty by definition (sampling floor). Also note `COALESCENCE_TIME` ($13.5$ samples) is shorter than the filter kernel ($17$ samples), so the two are entangled and cannot be tuned independently.

### D6. Terminology — "Omori law" — `OPEN`
$P(\Delta t) \sim \Delta t^{-1}$ is a claim about the **gap distribution**. The Omori law is a claim about the **rate**. These coincide only if the rate decays as $1/t$. Whether it does here is unresolved — the raw per-window counts rise steeply, but C1/D1 show those counts cannot be read as a rate without an exposure correction. **Until the rate is computed properly with exposure, do not call the result "Omori" in `project_summary.md` or anywhere else.**

### D7. Ringing frequency will change in the new rig — `OPEN`
If the ringing is a string mode, $f \propto 1/L$, so a longer pulley changes it and `RINGING_FREQ_HZ` needs re-measuring before any data is collected on the new setup.

---

## E. Minor conventions

1. Reference lines on the waiting-time PDF: $t^{-1}$ anchored to the histogram at $x_{\text{ref}} = 2$ s; $t^{-2}$ tail anchored to join it continuously at $t = 20$ s
2. `go.Scatter`, not `Scattergl` (WebGL unsupported on this machine)
3. All-windows overlay: rolling mean in log space, window 5 bins, `center=True`, `min_periods=1` — **display only, never fitted**
4. $P(M)$ drawn as dots, not bars (discrete data on a log axis)
5. $N = 10$ reliability line on the avalanche-count bar chart
6. `event_data.csv` is the deliverable to the supervisor — post-detection, one row per event, raw linear values
7. Aging-scatter filter `results_df["time"] > np.exp(3.5)`, scatter panel only; later set to `np.exp(0)`, i.e. effectively disabled
8. Bootstrap bands and Shohat reference lines on $P(M)$: **proposed and rejected**
9. Plotly subplot annotations: first subplot uses `"x domain"` / `"y domain"` unnumbered; `"x1 domain"` is rejected by Plotly

---

## F. Headline results — corrected status

| Quantity | Value | Status |
|---|---|---|
| $P(\Delta t)$ | $\sim \Delta t^{-1}$ | Holds. **Do not call it Omori** (D6) |
| Aging in shape | five windows, $\ln t = 2$–$7$: $-0.95, -1.00, -1.18, -1.19, -1.11$ | **Not a result** — no uncertainty (D4), censoring bias unquantified (D1) |
| Cutoff growth | $\max\Delta t \sim t^{0.98}$ | **Suspected artifact** (D3) — untested |
| Threshold type | fixed, 10 s | **Working default, not decided** (B1) |
| $P(M)$ | $\sim M^{-2.07}$ | **Parameter-dependent** — report $\alpha_M(\Delta t_{\text{thr}})$, not a number (B1) |
