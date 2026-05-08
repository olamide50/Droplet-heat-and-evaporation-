"""
plot_results.py – Visualise the droplet heating and evaporation results.

Generates 4 figures saved as PNG files:
  fig1_d2_law.png         – D²-law: normalised diameter² vs time (Cases 1 & 2)
  fig2_temperature.png    – Droplet temperature vs time (Cases 1 & 2)
  fig3_evap_rates.png     – Evaporation rate and heat flux vs time (Case 1)
  fig4_parametric.png     – Parametric sweep: K_D2 and lifetime vs T_inf

Usage:
    python3 plot_results.py
"""

import csv
import math
import os
import sys
import matplotlib
matplotlib.use("Agg")          # headless – no display required
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "lines.linewidth": 2.0,
    "figure.dpi": 130,
})

# ─────────────────────────────────────────────────────────────────────────────
# Helper: load a CSV into a dict of lists
# ─────────────────────────────────────────────────────────────────────────────
def load_csv(path):
    if not os.path.exists(path):
        sys.exit(f"ERROR: {path} not found – run ./droplet_sim first")
    rows = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            for k, v in row.items():
                rows.setdefault(k, []).append(float(v))
    return rows

c1 = load_csv("results_case1.csv")
c2 = load_csv("results_case2.csv")
pm = load_csv("results_parametric.csv")

# Convenience: time in ms, diameter in µm
def ms(col): return [v * 1e3 for v in col]
def us(col): return col   # already in µm in the CSV

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 – D²-law  (d/d0)²  vs  time
# ─────────────────────────────────────────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(7, 5))
ax1.plot(ms(c1["time_s"]), c1["d2_norm"], "b-",  label="Stagnant (U = 0 m/s)")
ax1.plot(ms(c2["time_s"]), c2["d2_norm"], "r--", label="Convective (U = 1 m/s)")

# Analytical D²-law line for Case 1 (linear fit through non-frozen zone)
t_arr = c1["time_s"]
d2_arr = c1["d2_norm"]
# Use last 60 % of points for the quasi-steady slope
n_qs = max(2, int(0.4 * len(t_arr)))
t_qs  = t_arr[n_qs:]
d2_qs = d2_arr[n_qs:]
if len(t_qs) >= 2:
    slope = (d2_qs[-1] - d2_qs[0]) / (t_qs[-1] - t_qs[0] + 1e-30)
    inter = d2_qs[0] - slope * t_qs[0]
    t_line = [t_arr[0], t_arr[-1]]
    d_line = [slope * t + inter for t in t_line]
    ax1.plot(ms(t_line), d_line, "b:", linewidth=1.2, label="D²-law fit (stagnant)")

ax1.set_xlabel("Time  [ms]")
ax1.set_ylabel(r"$(d/d_0)^2$")
ax1.set_title("D²-Law: Normalised Diameter Squared vs Time\n"
              "n-Heptane, d₀ = 100 µm, T∞ = 1000 K, P = 1 atm")
ax1.legend()
ax1.set_xlim(left=0)
ax1.set_ylim(0, 1.05)
ax1.grid(True, alpha=0.35)
fig1.tight_layout()
fig1.savefig("fig1_d2_law.png")
print("Saved fig1_d2_law.png")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 – Temperature evolution
# ─────────────────────────────────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(7, 5))
ax2.plot(ms(c1["time_s"]), c1["temperature_K"], "b-",  label="Stagnant (U = 0 m/s)")
ax2.plot(ms(c2["time_s"]), c2["temperature_K"], "r--", label="Convective (U = 1 m/s)")

# Wet-bulb lines
T_wb1 = c1["temperature_K"][-1]
T_wb2 = c2["temperature_K"][-1]
t_max = max(max(c1["time_s"]), max(c2["time_s"])) * 1e3
ax2.axhline(T_wb1, color="b", linestyle=":", linewidth=1.2,
            label=f"T_wb stagnant = {T_wb1:.1f} K")
ax2.axhline(T_wb2, color="r", linestyle=":", linewidth=1.2,
            label=f"T_wb convective = {T_wb2:.1f} K")
ax2.axhline(300.0,  color="gray", linestyle="-.", linewidth=1.0,
            label="T_d0 = 300 K")

ax2.set_xlabel("Time  [ms]")
ax2.set_ylabel("Droplet temperature  [K]")
ax2.set_title("Droplet Temperature vs Time\n"
              "n-Heptane, d₀ = 100 µm, T∞ = 1000 K, P = 1 atm")
ax2.legend(fontsize=9)
ax2.set_xlim(left=0)
ax2.grid(True, alpha=0.35)
fig2.tight_layout()
fig2.savefig("fig2_temperature.png")
print("Saved fig2_temperature.png")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 – Evaporation rate and heat flux (Case 1)
# ─────────────────────────────────────────────────────────────────────────────
fig3, (ax3a, ax3b) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)

mdot_ug_s = [abs(v) * 1e9 for v in c1["mdot_kg_s"]]   # µg/s
ax3a.plot(ms(c1["time_s"]), mdot_ug_s, "b-")
ax3a.set_ylabel("|ṁ|  [µg/s]")
ax3a.set_title("Evaporation Rate and Heat Transfer (Stagnant, T∞ = 1000 K)")
ax3a.grid(True, alpha=0.35)

Q_mW = [v * 1e3 for v in c1["Q_conv_W"]]               # mW
ax3b.plot(ms(c1["time_s"]), Q_mW, "r-", label="Q_conv")
ax3b.axhline(0, color="gray", linewidth=0.8)
ax3b.set_xlabel("Time  [ms]")
ax3b.set_ylabel("Q_conv  [mW]")
ax3b.grid(True, alpha=0.35)

ax3a.set_xlim(left=0)
fig3.tight_layout()
fig3.savefig("fig3_evap_rates.png")
print("Saved fig3_evap_rates.png")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 – Parametric sweep: K_D2 and lifetime vs T_inf
# ─────────────────────────────────────────────────────────────────────────────
fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(11, 5))

T_inf_list = pm["T_inf_K"]
K_list     = [v * 1e6 for v in pm["K_D2_m2s"]]   # mm²/s
t_list     = [v * 1e3 for v in pm["t_evap_s"]]    # ms

ax4a.plot(T_inf_list, K_list, "bs-", markersize=6)
ax4a.set_xlabel("Ambient temperature T∞  [K]")
ax4a.set_ylabel("Evaporation constant K  [mm²/s]")
ax4a.set_title("D²-Law Constant vs Ambient Temperature\n(n-Heptane, stagnant, d₀ = 100 µm)")
ax4a.grid(True, alpha=0.35)

ax4b.plot(T_inf_list, t_list, "rs-", markersize=6)
ax4b.set_xlabel("Ambient temperature T∞  [K]")
ax4b.set_ylabel("Droplet lifetime  [ms]")
ax4b.set_title("Droplet Lifetime vs Ambient Temperature\n(n-Heptane, stagnant, d₀ = 100 µm)")
ax4b.grid(True, alpha=0.35)

fig4.tight_layout()
fig4.savefig("fig4_parametric.png")
print("Saved fig4_parametric.png")

# ─────────────────────────────────────────────────────────────────────────────
print("\nAll plots saved.  Open the .png files to view the results.")
