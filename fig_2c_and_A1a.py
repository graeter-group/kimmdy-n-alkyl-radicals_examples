"""Script for analysis and plotting of figure 2c and A1a"""

import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
import experimental
import kimmdy_paper_theme
from kimmdy_paper_theme import single_column
from pathlib import Path
from matplotlib.transforms import Bbox
import util
plot_colors = kimmdy_paper_theme.auto_init()

cwd = Path(__file__).parent
T = 500

np.random.seed(0)

output_dir_root = cwd / "output" / "n-alkyl_examples"

rates = np.load(output_dir_root / 'rates_all.npy', allow_pickle='TRUE').item()
rates = {shift: {alkyl: rs["r"] for alkyl, rs in shift_dict.items()} for shift, shift_dict in rates.items() }
shifts = sorted(rates.keys())

identity_shifts = {
    "1-2": "ethyl", "1-3": "propyl", "1-4": "butyl", "1-5": "pentyl", "1-6": "hexyl", "1-7": "heptyl", "1-8": "octyl",
}

alkyls = list(rates["1-2"].keys())
rates_by_alkyl = {alkyl: {} for alkyl in alkyls}
for alkyl in alkyls:
    for shift, shift_rates in rates.items():
        if len(rates[shift][alkyl]) == 0: continue
        tmp = np.array([t for t in rates[shift][alkyl]])
        rates_by_alkyl[alkyl].update({shift: tmp})
pprint(rates_by_alkyl)

relative_rate_alkyls = ["heptyl", "octyl"]
relative_rate_shifts = ["1-2", "1-3", "1-4", "1-5", "1-6", "1-7",]

relative_rates_by_alkyl = {alkyl: {shift: rates_by_alkyl[alkyl][shift] for shift in relative_rate_shifts} for alkyl in relative_rate_alkyls}
R = np.zeros((len(relative_rate_alkyls), len(relative_rate_shifts), 20))
for i, (alkyl, alkyl_dict) in enumerate(relative_rates_by_alkyl.items()):
    for j, (shift, r) in enumerate(alkyl_dict.items()):
        if (shift == "1-2") or (shift == "1-7"): r = 0
        R[i, j, :] = r

# Rates divided by sum of rates over all shifts -> for each run new probability (20)
P = R / np.sum(R, axis=1)[:, np.newaxis, :]

assert np.allclose(np.sum(P, axis=1), np.ones((R.shape[0], R.shape[-1])))

mean_rates_per_shift = {}
for shift in shifts:
    kmdy_shift_rates = rates[shift]
    rs_alkyl_mean_rate = {}
    for alkyl, rs in kmdy_shift_rates.items():
        if alkyl in identity_shifts[shift] or len(rs) == 0: continue
        rs_alkyl_mean_rate[alkyl] = np.mean(rs)
    if len(rs_alkyl_mean_rate) > 0:
        mean_rates_per_shift[shift] = rs_alkyl_mean_rate

pprint(mean_rates_per_shift)
mean_rates_per_over_all_molecules = {shift: float(np.mean([r for alkyl, r in alkyl_rates.items()])) for shift, alkyl_rates in
                                     mean_rates_per_shift.items()}
indiv_rates_per_over_all_molecules = {shift: np.array([r for alkyl, r in alkyl_rates.items()]) for shift, alkyl_rates in
                                      mean_rates_per_shift.items()}

pprint(mean_rates_per_over_all_molecules)
pprint(indiv_rates_per_over_all_molecules)

fig, axes = plt.subplots(2, 1, figsize=(single_column, single_column), sharex=False)
width = 1 / 4
for i, ax in enumerate(axes):

    shifts = sorted(mean_rates_per_over_all_molecules.keys())
    xticks = {shift: i for i, shift in enumerate(shifts)}

    exp_rates = {
        shift: np.array([validation_paper_result.get_rate(T) for validation_paper_result in experimental.rate_results_all[shift]])
        for shift in shifts if shift not in ["1-2", "1-7", "1-8"]
    }
    exp_rates = exp_rates | {"1-2": [0.0], "1-7": [0.0]}
    pprint(exp_rates)

    exp_shifts = sorted(exp_rates.keys())
    exp_shift_means = np.array([np.mean(exp_rates[shift]) for shift in exp_shifts])
    exp_shift_values = [exp_rates[shift] for shift in exp_shifts]

    if i == 0:
        xticks_relative = {shift: x for shift, x in xticks.items() if shift in relative_rate_shifts}

        p_experimental = exp_shift_means / sum(exp_shift_means)
        ax.bar(np.array([xticks[shift] + 3/2 * width for shift in exp_shifts]), p_experimental,
               alpha=1.0, width=width,
               color=plot_colors["experiment"], zorder=1, label="Experimental: n-alkyl")

        colors = {"heptyl": plot_colors["kimmdy_light"], "octyl": plot_colors["kimmdy"]}

        relabel = {"heptyl": "Heptyl", "octyl": "Octyl"}

        for i_alkyl, relative_alkyl in enumerate(relative_rate_alkyls):
            for j_shift, shift in enumerate(relative_rate_shifts):
                if shift in ["1-2", "1-7"]: continue
                ps_alkyl = P[i_alkyl, j_shift]
                ax.scatter([xticks[shift]+(i_alkyl-1/2) * width]*ps_alkyl.size+ np.random.normal(0, 0.01, len(ps_alkyl)), ps_alkyl,
                           alpha=0.8,
                           color="white", zorder=100, facecolor="None", edgecolors="black", linewidth=0.5)

            p_alkyl = np.mean(P[i_alkyl], axis=1)
            brier_kimmdy = np.mean([util.brier_score_divergence(p, p_experimental) for p in P[i_alkyl].T])
            brier_ignorance = util.brier_score_divergence(np.ones(p_experimental.size)/p_experimental.size, p_experimental)

            print(f"{relative_alkyl} - BD: {brier_kimmdy:.2f}, IG: {brier_ignorance:.2f}")

            ax.bar(np.array(list(xticks_relative.values())) + (i_alkyl-1/2) * width, p_alkyl, alpha=1.0, width=width,
                   zorder=1,
                   label=f"KIMMDY: {relabel[relative_alkyl]}",
                   color=colors[relative_alkyl],)

        ax.set_ylim((0, 1))
        ax.set_ylabel("Selection Probability")

        ax.set_xticks(np.arange(len(shifts)) + 0.125, shifts)
        ax.set_xlabel("Hydrogen Atom Transfer")

    else:
        shift_means = [mean_rates_per_over_all_molecules[shift] for shift in shifts]
        shift_values = [indiv_rates_per_over_all_molecules[shift] for shift in shifts]

        ax.bar(xticks.values(), shift_means, alpha=1.0, width=width, color=plot_colors["kimmdy"], zorder=1, label="KIMMDY: n-alkyl")

        for j, (shift, shift_vs) in enumerate(zip(shifts, shift_values)):
            ax.scatter([xticks[shift]] * len(shift_vs) + np.random.normal(0, 0.01, len(shift_vs)), shift_vs, alpha=0.8,
                       color="white", zorder=1, facecolor="None", edgecolors="black", linewidth=0.5)

        ax.bar(np.array([xticks[shift] for shift in exp_shifts]) + width, exp_shift_means, alpha=1.0, width=width,
               color=plot_colors["experiment"], zorder=1, label="Experimental: n-alkyl")

        for j, (exp_shift, exp_shift_vs) in enumerate(zip(exp_shifts, exp_shift_values)):
            ax.scatter([xticks[exp_shift] + width] * len(exp_shift_vs) + np.random.normal(0, 0.01, len(exp_shift_vs)),
                       exp_shift_vs, alpha=0.5,
                       color="white", zorder=1, edgecolors="black", linewidth=0.5)

        ax.set_xticks(np.arange(len(shifts)) + 0.125, shifts)
        ax.set_ylabel(f"Absolute rates at {T}K [1/s]")
        ax.set_xlabel("Hydrogen Atom Transfer")

        if i == 1:
            ax.set_yscale("log")
        else:
            legend = ax.legend(loc="upper left", framealpha=0)

        ax.legend(loc='upper left', frameon=False, framealpha=0)

axes[1].set_xlim(axes[1].get_xlim()[0], xticks["1-7"]+width+0.45)
axes[0].legend(loc='upper left', frameon=False, framealpha=0)

# Upper figure
plt.tight_layout(w_pad=0.0, h_pad=0.0, pad=0.5)
default_extent = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
x0, y0, x1, y1 = default_extent.x0, default_extent.y0, default_extent.x1, default_extent.y1

# left, bottom, right, top
custom_bbox = Bbox.from_extents(x0+0.1, y1/2, x1-0.05, y1-0.05)
plt.savefig(output_dir_root / "rates_relative.png", bbox_inches=custom_bbox, dpi=300)
kimmdy_paper_theme.convert_to_rgb(output_dir_root / "rates_relative.png")


# Lower figure
custom_bbox = Bbox.from_extents(x0+0.1, y0+0.1, x1-0.05, y1/2)
plt.savefig(output_dir_root / "rates_absolute.png", bbox_inches=custom_bbox, dpi=300)
kimmdy_paper_theme.convert_to_rgb(output_dir_root / "rates_absolute.png")

plt.show()
