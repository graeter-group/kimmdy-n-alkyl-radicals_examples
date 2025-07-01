"""Script for extracting reaction rates from n-alkyl simulation runs (basis for fig 2c and A1a)"""

from pathlib import Path
import numpy as np
import util
from pprint import pprint, pformat

cwd = Path(__file__).parent

main_run_folder = cwd / "simulations" / "n-alkyl_examples" / "systems" / "500K"

stride = 100
T = 500
start_t = 0
R = 1.987e-3
h_planck = 4.135667696e-15

output_dir_root = cwd / "output" / "n-alkyl_examples"
output_dir_root.mkdir(exist_ok=True)

radical_indices = {
    "propyl": dict(
        C_rad=3,
        shift_indices={
            "1-2": (7, 8),
            "1-3": (4, 5, 6),
        }
    ),
    "butyl": dict(
        C_rad=4,
        shift_indices={
            "1-2": (10, 11),
            "1-3": (8, 9),
            "1-4": (5, 6, 7),
        }
    ),
    "pentyl": dict(
        C_rad=5,
        shift_indices={
            "1-2": (13, 14),
            "1-3": (11, 12),
            "1-4": (9, 10),
            "1-5": (6, 7, 8),
        }
    ),
    "hexyl": dict(
        C_rad=6,
        shift_indices={
            "1-2": (16, 17),
            "1-3": (14, 15),
            "1-4": (12, 13),
            "1-5": (10, 11),
            "1-6": (7, 8, 9),
        }
    ),
    "heptyl": dict(
        C_rad=7,
        shift_indices={
            "1-2": (19, 20),
            "1-3": (17, 18),
            "1-4": (15, 16),
            "1-5": (13, 14),
            "1-6": (11, 12),
            "1-7": (8, 9, 10),
        }
    ),
    "octyl": dict(
        C_rad=8,
        shift_indices={
            "1-2": (22, 23),
            "1-3": (20, 21),
            "1-4": (18, 19),
            "1-5": (16, 17),
            "1-6": (14, 15),
            "1-7": (12, 13),
            "1-8": (9, 10, 11),
        }
    ),
}

molecules = list(radical_indices.keys())


order_main = ["1-2", "1-3", "1-4", "1-5", "1-6", "1-7", "1-8"]
all_mean_rates = {shift: [] for shift in order_main}

rates_all = {shift: {m: {"r": []} for m in molecules} for shift in order_main}  # shift, molecule, index


for alkyl in radical_indices.keys():
    alkyl_run_folders = sorted([f.parent for f in main_run_folder.glob(f"{alkyl}/*/kimmdy.history")])
    output_dir = output_dir_root / alkyl
    output_dir.mkdir(parents=True, exist_ok=True)

    n_runs = len(alkyl_run_folders)

    assert n_runs == 20, f"{alkyl} has {n_runs} runs"

    order = [o for o in order_main if o in radical_indices[alkyl]["shift_indices"].keys()]

    rates_kimmdy_dict = {o: [] for o in order}

    for i_run, run_folder in enumerate(alkyl_run_folders):

        run_name = run_folder.name
        run_output_dir = output_dir / run_folder.name
        run_output_dir.mkdir(parents=True, exist_ok=True)

        rates_csv = run_folder / "3_decide_recipe" / "recipes.csv"
        se_dir = run_folder / "2_hat_reaction" / "se"

        df = util.load_results(rates_csv, se_dir, T, stride=stride, start_t=start_t)

        C_rad = radical_indices[alkyl]["C_rad"]
        shift_indices = radical_indices[alkyl]["shift_indices"]

        inv_shift_indices = {i: k for k, v in shift_indices.items() for i in v}
        df["ids"] = df["ids"].apply(lambda x: [int(x_) for x_ in x])
        assert df["ids"].apply(lambda x: C_rad in x).all()
        assert df["ids"].apply(lambda x: len(x) == 2).all()
        df["index_H"] = df["ids"].apply(lambda x: x[0] if x[1] == C_rad else x[1])
        df["shift"] = df["index_H"].apply(lambda x: inv_shift_indices[x])

        df.to_csv(run_output_dir / f"rates.csv")

        t_tot = df["time (s)"].max() - df["time (s)"].min()

        rates_per_index = {idx:{"r_kmdy": np.nan, "shift": inv_shift_indices[idx]}
                           for idx in sorted(df["index_H"].unique())}
        run_rates = {shift: [] for shift in df["shift"].unique()}
        for index_H, index_H_group in df.groupby("index_H"):
            shift_H = inv_shift_indices[index_H]
            mean_rate_H = index_H_group["rates (1/s)"].sum() * stride * 1e-15 / t_tot

            rates_per_index[index_H]["r_kmdy"] = mean_rate_H

            n_shift_H = len(shift_indices[shift_H])
            run_rates[shift_H].append(mean_rate_H)

        for shift in df["shift"].unique():
            rates_all[shift][alkyl]["r"].append(np.sum(run_rates[shift]))

        rates_kimmdy_run = {}
        for i_shift, (shift, shift_group) in enumerate(df.groupby("shift")):
            mean_rate = shift_group["rates (1/s)"].sum() * stride * 1e-15 / t_tot

            rates_kimmdy_dict[shift].append(mean_rate)
            print(f"{alkyl}: Shift: {shift}; Mean rate: {mean_rate:.1E}")

    rates_kimmdy = [[rates_kimmdy_dict[o][i_run] for o in order] for i_run in range(n_runs)]
    rates_kimmdy = np.array(rates_kimmdy)

    p_rates_kimmdy = (rates_kimmdy.T / np.sum(rates_kimmdy, axis=1)).T

    mean_rates = np.mean(rates_kimmdy, axis=0)

    res = {shift: f"{mean_rates[i]:0.1e}" for i, shift in enumerate(shift_indices.keys())}
    with open(output_dir / 'rates.txt', 'w') as f:
        f.write(pformat(res, indent=1))

    for i, o in enumerate(order):
        all_mean_rates[o].append(mean_rates[i])

res_all = {shift: [f"{x:0.1e}" for x in all_mean_rates[shift]] for i, shift in enumerate(order_main)}
with open(output_dir_root / 'all_rates.txt', 'w') as f:
    f.write(pformat(res_all, indent=1))


pprint(rates_all)
with open(output_dir_root / 'all_rates_individual.txt', 'w') as f:
    f.write(pformat(rates_all, indent=1))
np.save(output_dir_root / 'rates_all.npy', rates_all)
