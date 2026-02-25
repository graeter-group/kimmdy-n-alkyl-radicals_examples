"""Script for analysis and plotting of figure 3d and table S5"""
from pathlib import Path
import numpy as np
import pandas as pd
import util
import matplotlib.pyplot as plt
from pprint import pprint
import experimental
import kimmdy_paper_theme
from kimmdy_paper_theme import single_column
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import Bbox
plot_colors = kimmdy_paper_theme.auto_init()

cwd = Path(__file__).parent

parameter_results = list((cwd / "output" / "hyperparameters" / "simulations" / "hyperparameter_grid" / "systems" / "octyl" / "parameter_runs").glob("n_steps-*__n_unique-*"))

parameter_rates = {}
t_sims = set()
n_uniques = set()
n_runs = -1
for p in parameter_results:
    rates = np.load(p / "rates_all.npy", allow_pickle=True).item()
    rates = {shift: np.array(list(alkyl_dict["octyl"].values())).squeeze().tolist() for shift, alkyl_dict in rates.items()
             if shift not in ["1-8"]}

    t_sim = int(p.stem.split("__")[0].split("n_steps-")[1])
    n_unique = int(p.stem.split("__")[1].split("n_unique-")[1])

    if n_unique == 100000: continue  # did not finish running

    if n_runs == -1:
        n_runs = len(rates["1-2"])

    for shift, rs_shift in rates.items():
        if type(rs_shift) is float:
            rs_shift = [rs_shift]
        if len(rs_shift) != n_runs:
            rs_new = [np.nan] * n_runs
            for i_run in range(len(rs_shift)):
                rs_new[i_run] = rs_shift[i_run]
            rates[shift] = rs_new

    if t_sim not in parameter_rates:
        parameter_rates[t_sim] = {}
    parameter_rates[t_sim][n_unique] = rates
    t_sims.add(t_sim)
    n_uniques.add(n_unique)

pprint(parameter_rates)

t_sims = np.array(sorted(list(t_sims))).astype(int)
n_uniques = np.array(sorted(list(n_uniques))).astype(int)

T = 500

shifts = sorted(rates.keys())

exp_rates = {
        shift: np.array([validation_paper_result.get_rate(T) for validation_paper_result in experimental.rate_results_all[shift]]).tolist()
        for shift in experimental.rate_results_all.keys()
    }
exp_shifts = exp_rates.keys()
pprint(exp_rates)

r_experimental = np.array([np.mean(exp_rates[shift]) for shift in exp_shifts])
p_experimental = r_experimental / r_experimental.sum()

R_longest_scatter = np.full((n_uniques.size, len(shifts), 4), -1.0)  # All rates for longest simlation time

BD = np.full((t_sims.size, n_uniques.size), np.nan)
R_mean = np.full((len(shifts), t_sims.size, n_uniques.size), -1.0)
R_min = np.full((len(shifts), t_sims.size, n_uniques.size), -1.0)
R_max = np.full((len(shifts), t_sims.size, n_uniques.size), -1.0)
assert len(shifts) == 6

ratios = {}
for i, t_sim in enumerate(t_sims):
    for j, n_unique in enumerate(n_uniques):

        if n_unique > t_sim: continue
        rs_kimmdy = parameter_rates[int(t_sim)][int(n_unique)]

        rs_kimmdy = np.array([rs_kimmdy[shift] for shift in shifts])
        rs_kimmdy[np.isnan(rs_kimmdy)] = 0.0

        ratio = float(n_unique / t_sim)
        if ratio not in ratios.keys():
            ratios[ratio] = []

        R_mean[:, i, j] = np.nanmean(rs_kimmdy, axis=1)
        R_min[:, i, j] = np.nanmin(rs_kimmdy, axis=1)
        R_max[:, i, j] = np.nanmax(rs_kimmdy, axis=1)

        # For BD analysis, we exclude 1-2 and 1-7
        assert rs_kimmdy.shape[0] == 6
        assert p_experimental.shape[0] == 4
        p_kimmdy = (rs_kimmdy[1:-1] / np.nansum(rs_kimmdy[1:-1], axis=0)).T

        if t_sim == 100000 and n_unique == 100000:
            print(t_sim, n_unique)
            print(repr(p_kimmdy))
            print(p_kimmdy.shape, shifts)

        if t_sim == 100000:
            R_longest_scatter[j, :, :] = rs_kimmdy

        bd = np.mean([util.brier_score_divergence(p_k, p_experimental) for p_k in p_kimmdy])
        assert np.isclose(p_experimental.sum(), 1.0) and np.isclose(p_kimmdy[0].sum(), 1.0)

        print(bd)

        BD[i, j] = bd

data = np.hstack([t_sims[..., None], BD])
df = pd.DataFrame(data, columns=np.hstack([r"n\textsubscript{frames}", n_uniques]))
print(df.map(lambda x: f"{x:.2f}").fillna('X').to_latex(index=False).replace("nan", "X"))

custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", ["#C3006B", "#019050"])

shift_colors = {shift: custom_cmap(i/len(shifts)) for i, shift, in enumerate(shifts)}

R_longest = R_mean[:, -1, :]  # shift, n_unique

N, M = R_longest.shape
x = np.arange(M)*1.1  # M groups -> # of predictions
width = 1 / N   # width of each bar

# Create plot
fig, ax = plt.subplots(figsize=(single_column, single_column/2))

df_scatter_rates_per_shift = {}

for i in range(N):  # iterate over shifts: first build 1-2, afterwards 1-3 etc
    bars = ax.bar(x + i*width, R_longest[i], width, label=shifts[i], color=shift_colors[shifts[i]])
    bar = bars[0]
    x_center = bar.get_x() + bar.get_width() / 2
    # Plot also the individual run rates
    ax.scatter(np.repeat(x + i*width, 4), R_longest_scatter[:, i, :].flatten(), alpha=0.5,
               color="white", zorder=1, facecolor="None", edgecolors="black", linewidth=0.5)
    ax.text(
        x_center+0.03,
        bar.get_height()+bar.get_height()*7,
        str(shifts[i]),
        ha='center',
        va='bottom',
        rotation=90,
        fontsize=6.5,
        color=shift_colors[shifts[i]],
        weight='bold',
    )

    # Save for later exporting
    df_shift = pd.DataFrame(
            data=R_longest_scatter[:, i, :],
            index=n_uniques,
            columns=[f"Rate run {k} [1/s]" for k in range(4)]
        )
    df_shift.index.name = "# of predicted reactions"
    df_scatter_rates_per_shift[shifts[i]] = df_shift
with pd.ExcelWriter(cwd / "output" / "hyperparameters" / f"alkyl_rate_convergence.xlsx") as writer:
    for shift, df in df_scatter_rates_per_shift.items():
        df.to_excel(writer, sheet_name=f"fig_3d_shift_{shift}", index=True)

# Labeling
ax.set_xticks(x + width*(N-1)/2)
ax.set_xticklabels(n_uniques)
ax.set_ylabel('Rate [1/s]')
ax.set_xlabel('Number of predicted reactions')
plt.tight_layout()

plt.yscale("log")

fig.subplots_adjust(top=1.0, bottom=0, left=0, right=1)

plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)

default_extent = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
x0, y0, x1, y1 = default_extent.x0, default_extent.y0, default_extent.x1, default_extent.y1

# left, bottom, right, top
custom_bbox = Bbox.from_extents(x0, y0, x1, y1)

plt.savefig(cwd / "output" / "hyperparameters" / "alkyl_rate_convergence.png", bbox_inches=custom_bbox, dpi=300)
plt.savefig(cwd / "output" / "hyperparameters" / "alkyl_rate_convergence.pdf", bbox_inches=custom_bbox, dpi=300)

plt.show()
