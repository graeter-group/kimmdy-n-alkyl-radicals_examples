"""Script for extracting reaction rates from simulation runs for different combinations of t_sim and n_unique (basis for fig 2d and A1b)"""

from pathlib import Path
import numpy as np
import time
import util

cwd = Path(__file__).parent

top_folder = cwd / "simulations" / "hyperparameter_grid" / "systems" / "octyl" / "parameter_runs"
subfolders = top_folder.glob("n_steps-*__n_unique-*")

for subfolder in subfolders:
    main_run_folder = subfolder
    run_alkyl = main_run_folder.parent.parent.name

    stride = 100
    T = 500
    start_t = 0

    output_dir_root = cwd / "output" / "hyperparameters" / main_run_folder.relative_to(cwd)

    radical_indices = {
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

    assert run_alkyl in radical_indices.keys()

    molecules = list(radical_indices.keys())

    order_main = ["1-2", "1-3", "1-4", "1-5", "1-6", "1-7", "1-8"]
    all_mean_rates = {shift: [] for shift in order_main}

    rates_all = {shift: {m: {"r": []} for m in molecules} for shift in order_main}  # shift, molecule, index

    for alkyl in radical_indices.keys():
        if alkyl != run_alkyl: continue
        alkyl_run_folders = sorted([f.parent for f in main_run_folder.glob(f"*/kimmdy.history")])
        output_dir = output_dir_root / alkyl
        output_dir.mkdir(parents=True, exist_ok=True)

        n_runs = len(alkyl_run_folders)
        if n_runs == 0:
            print(f"{alkyl} has {n_runs} runs. SKIPPING!")
            time.sleep(1)
            continue

        assert n_runs == 4, f"{alkyl} has {n_runs} runs"

        order = [o for o in order_main if o in radical_indices[alkyl]["shift_indices"].keys()]

        rates_kimmdy_dict = {o: [np.nan]*n_runs for o in order}
        found_rates = False

        for i_run, run_folder in enumerate(alkyl_run_folders):

            run_name = run_folder.name
            run_output_dir = output_dir / run_folder.name
            run_output_dir.mkdir(parents=True, exist_ok=True)

            rates_csv = run_folder / "2_decide_recipe" / "recipes.csv"
            se_dir = run_folder / "1_hat_reaction" / "se"
            if not rates_csv.exists():
                print(f"{rates_csv} does not exist")
                continue

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

            rates_per_index = {idx: {"r_kmdy": np.nan, "shift": inv_shift_indices[idx]}
                               for idx in sorted(df["index_H"].unique())}
            run_rates = {shift: [] for shift in df["shift"].unique()}
            for i_index_H, (index_H, index_H_group) in enumerate(df.groupby("index_H")):
                shift_H = inv_shift_indices[index_H]
                n = int(subfolder.name.split('__')[0].split('-')[1])
                mean_rate_H = index_H_group["rates (1/s)"].sum() / n
                rates_per_index[index_H]["r_kmdy"] = mean_rate_H

                run_rates[shift_H].append(mean_rate_H)

            for shift in df["shift"].unique():
                rates_all[shift][alkyl]["r"].append(np.sum(run_rates[shift]))

            rates_kimmdy_run = {}
            for i_shift, (shift, shift_group) in enumerate(df.groupby("shift")):
                n = int(subfolder.name.split('__')[0].split('-')[1])
                mean_rate = shift_group["rates (1/s)"].sum() / n
                rates_kimmdy_dict[shift][i_run] = float(mean_rate)

            found_rates = True

        if not found_rates:
            continue

        rates_kimmdy = [[rates_kimmdy_dict[o][i_run] for o in order] for i_run in range(n_runs)]
        rates_kimmdy = np.array(rates_kimmdy)

        p_rates_kimmdy = (rates_kimmdy.T / np.nansum(rates_kimmdy, axis=1)).T

        mean_rates = np.mean(rates_kimmdy, axis=0)

        res = {shift: f"{mean_rates[i]:0.1e}" for i, shift in enumerate(shift_indices.keys())}

        for i, o in enumerate(order):
            all_mean_rates[o].append(mean_rates[i])

    np.save(output_dir_root / 'rates_all.npy', rates_all)