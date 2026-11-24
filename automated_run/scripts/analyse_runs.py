#%%
from pathlib import Path
import os
import ast
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from tqdm import tqdm

# Override with ANALYSIS_ROOT when running on the cluster (local sshfs path is the default)
root_dir = Path(os.environ.get(
    "ANALYSIS_ROOT",
    "/home/hartmanne/Cluster/otter_ptmp/alkyl_simulations/kimmdy-n-alkyl-radicals_examples/automated_run/",
))
script_dir = root_dir / "scripts"
# Force this directory to the top of Python's search path
if script_dir.as_posix() not in sys.path:
    sys.path.insert(0, script_dir.as_posix())

# Import rate results from experimental.py
from experimental import rate_results, R_cal, R_kcal

#%%
# --- Global Configurations ---
TARGET_TEMP_K = 500.0  # Temperature to evaluate experimental rate equations
output_dir = root_dir / "pngs"
#output_dir.mkdir(exist_ok=True)

# Override with ANALYSIS_RUN_NAME to analyse another run without editing the script
run_name = os.environ.get("ANALYSIS_RUN_NAME", "HAT-MACE_datasetv2_capdist8_20261122_085918")
run_dir = root_dir / "runs" / run_name

systems = ["heptyl", "octyl"]

def trajectory_frames(mdp_path: Path) -> tuple[int, int]:
    """Number of trajectory frames and frame spacing (fs) of the equilibrium MD.

    KIMMDY scans the .trr, so frames are written every nstxout steps.
    Frames not evaluated by KIMMDY (e.g. beyond n_unique per reaction) count as rate 0,
    so time averages must divide by the total number of trajectory frames.
    """
    params = {}
    for line in mdp_path.read_text().splitlines():
        line = line.split(";", 1)[0]
        if "=" in line:
            key, value = line.split("=", 1)
            params[key.strip().replace("_", "-")] = value.strip()
    nsteps, nstxout, dt_ps = int(params["nsteps"]), int(params["nstxout"]), float(params["dt"])
    return nsteps // nstxout + 1, int(round(nstxout * dt_ps * 1000))


file_paths: list[Path] = []
traj_info: dict[Path, tuple[int, int]] = {}  # recipes.csv -> (n_total_frames, frame_dt_fs)

for system in systems:
    sys_frames = trajectory_frames(run_dir / system / "sd.mdp")
    print(f"{system}: {sys_frames[0]} trajectory frames every {sys_frames[1]} fs")
    for fp in sorted((run_dir / system).glob("sim*/HAT_000/2_decide_recipe/recipes.csv")):
        file_paths.append(fp)
        traj_info[fp] = sys_frames

if len(file_paths) == 0:
    raise FileNotFoundError("No valid recipes.csv files were found in specified paths.")
#%% functions
# --- Helper Parsing & Data Processing Functions ---
def extract_ordered_ids(recipe_str: str) -> tuple[int, ...]:
    """Extracts unique atom index sequence from recipe strings."""
    raw_indices = re.findall(r"atom_ix_\d+=(\d+)", str(recipe_str))
    seen = set()
    unique_indices = []
    for idx_str in raw_indices:
        idx = int(idx_str)
        if idx not in seen:
            seen.add(idx)
            unique_indices.append(idx)
    return tuple(unique_indices)


def extract_first_timespan(timespan_str: str) -> float:
    """Extracts the first timespan value (ps). Kept as float: int() would truncate to whole ps."""
    val = ast.literal_eval(str(timespan_str))
    return float(val[0][0])


def extract_first_rate(rate_str: str) -> float:
    """Extracts the first rate constant value."""
    val = ast.literal_eval(str(rate_str))
    return float(val[0])


def load_and_parse(file_path: str) -> pd.DataFrame:
    """Loads CSV and converts raw strings into structured columns."""
    raw_df = pd.read_csv(file_path)
    return pd.DataFrame(
        {
            "ids": raw_df["recipe_steps"].apply(extract_ordered_ids),
            "rates": 1e12 * raw_df["rates"].apply(extract_first_rate),  # s⁻¹
            # ps -> fs, rounded to whole fs (not truncated to whole ps)
            "timespans": (1000 * raw_df["timespans"].apply(extract_first_timespan)).round().astype(np.int64),
        }
    )


def process_subsampling(
    df: pd.DataFrame, file_id: int, n_total_frames: int, frame_dt_fs: int
) -> pd.DataFrame:
    """Processes an already parsed simulation DataFrame for subsampling down to different frame resolutions.

    Subsampling keeps every step-th trajectory frame (by frame index, so frames without
    any evaluated recipe still count) and averages over all frames of the subsample.
    """
    total_sim_time_fs = (n_total_frames - 1) * frame_dt_fs
    total_sim_time_s = 1e-15 * total_sim_time_fs

    # Filter out valid HAT recipes (at least 3 atom IDs)
    df_valid = df[df["ids"].apply(lambda x: len(x) >= 3)].copy()
    df_valid["id_diff"] = df_valid["ids"].apply(lambda x: abs(x[0] - x[2]))

    # Trajectory frame index of each recipe
    frame_idx = (df_valid["timespans"] / frame_dt_fs).round().astype(np.int64)

    all_subsamples = []
    step_sizes = [1, 10, 100, 1000, 10000, 100000]

    for step in step_sizes:
        # Keep trajectory frames 0, step, 2*step, ...
        subsampled_df = df_valid[frame_idx % step == 0]

        n_frames = (n_total_frames - 1) // step + 1  # Number of trajectory frames in the subsample
        if subsampled_df.empty or n_frames == 0:
            continue

        print(f"Step: {step}, Selected Frames: {n_frames}")
        label = "1 (Full)" if step == 1 else f"1/{step}"

        # --- STEP 1: First grouping by unique 'ids' ---
        grouped_ids = (
            subsampled_df.groupby(["ids", "id_diff"], as_index=False)
            .agg(sum_rate=("rates", "sum"))
        )

        # True time-averaged rate constant over the sampled frames
        grouped_ids["rate"] = grouped_ids["sum_rate"] / n_frames  # s⁻¹

        # Reaction probability over trajectory: P ≈ <k> * total_sim_time_s
        grouped_ids["reaction_probability"] = grouped_ids["rate"] * total_sim_time_s

        # --- STEP 2: Second grouping to aggregate from 'ids' up to 'id_diff' ---
        grouped = (
            grouped_ids.groupby("id_diff", as_index=False)
            .agg(
                rate=("rate", "sum"),
                reaction_probability=("reaction_probability", "sum")
            )
        )

        # Normalize reaction probabilities to create selection_probability
        total_rxn_prob = grouped["reaction_probability"].sum()
        grouped["selection_probability"] = (
            grouped["reaction_probability"] / total_rxn_prob if total_rxn_prob > 0 else 0.0
        )

        # Metadata & reaction type labels
        grouped["file_id"] = file_id
        grouped["sampling_fraction"] = label
        grouped["step_size"] = step
        grouped["num_frames"] = n_frames
        grouped["reaction_type"] = grouped["id_diff"].apply(
            lambda d: f"1-{int(d) + 1}"
        )

        all_subsamples.append(grouped)

    if not all_subsamples:
        return pd.DataFrame()

    return pd.concat(all_subsamples, ignore_index=True)


# --- Experimental Data Extractor ---
def get_experimental_dataframe(temp_k: float = 500.0) -> pd.DataFrame:
    """Evaluates experimental RateResult records at a target temperature."""
    exp_records = []

    # Inject T and gas constants into global scope so eval() inside get_rate() works
    globals()["T"] = temp_k
    globals()["R_cal"] = R_cal
    globals()["R_kcal"] = R_kcal

    # Iterate over the INSTANCE LIST rate_results (not the class RateResult)
    for r in rate_results:
        try:
            k_val = r.get_rate(temp_k)
            exp_records.append(
                {
                    "reaction_type": r.reaction_type,
                    "base_molecule": r.base_molecule,
                    "paper_id": r.paper_id,
                    "rate": float(k_val),
                    "source": f"Exp (Paper {r.paper_id})",
                    "origin": "Experiment",
                }
            )
        except Exception as e:
            print(f"Skipping paper ID {r.paper_id} ({r.reaction_type}): {e}")

    df_exp = pd.DataFrame(exp_records)

    # Calculate selection probabilities grouped by paper source (paper_id)
    exp_probs = []
    for paper_id, g in df_exp.groupby("paper_id"):
        g_copy = g.copy()
        tot_rate = g_copy["rate"].sum()
        g_copy["reaction_probability"] = (
            g_copy["rate"] / tot_rate if tot_rate > 0 else 0.0
        )
        exp_probs.append(g_copy)

    return pd.concat(exp_probs, ignore_index=True) if exp_probs else df_exp

#%%
# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

print("Loading raw simulation files...")
# Load and parse raw files ONCE to store in memory
raw_parsed_dfs = []
full_dfs = []

for idx, fp in tqdm(enumerate(file_paths), desc="Parsing simulation files"):
    df_raw = load_and_parse(fp.as_posix())
    raw_parsed_dfs.append(df_raw)
    
    # 1. Get total frame count of the trajectory (including frames with 0 predicted rate
    #    and frames KIMMDY did not evaluate); len(df_raw) would count recipe rows instead
    n_total_frames, frame_dt_fs = traj_info[fp]
    total_sim_time_fs = (n_total_frames - 1) * frame_dt_fs
    total_sim_time_s = 1e-15*total_sim_time_fs
    # Recipe times are float32 (from the .trr): ~1e5 ps carries a few fs of rounding error,
    # so only flag times beyond the trajectory or clearly off the frame grid (>10% of spacing)
    frame_idx_all = (df_raw["timespans"] / frame_dt_fs).round()
    grid_offset_fs = (df_raw["timespans"] - frame_idx_all * frame_dt_fs).abs().max()
    if frame_idx_all.max() >= n_total_frames or grid_offset_fs > 0.1 * frame_dt_fs:
        print(f"WARNING {fp}: recipe times do not match {n_total_frames} frames every {frame_dt_fs} fs")
    
    # 2. Filter to valid HAT recipes
    df_valid = df_raw[df_raw["ids"].apply(lambda x: len(x) >= 3)].copy()
    df_valid["id_diff"] = df_valid["ids"].apply(lambda x: abs(x[0] - x[2]))
    
    # 3. Sum rates for active frames per unique ID
    grouped_ids = df_valid.groupby(["ids", "id_diff"], as_index=False).agg(
        sum_rate=("rates", "sum")
    )
    
    # 4. Average over ALL trajectory frames (active + zero frames)
    grouped_ids["rate"] = grouped_ids["sum_rate"] / n_total_frames  # s⁻¹
    
    # 5. Aggregate up to id_diff (reaction type)
    grouped = grouped_ids.groupby("id_diff", as_index=False).agg(
        rate=("rate", "sum")  # sum average rates of equivalent atom paths
    )
    
    # 6. Calculate total integrated probability for selection
    grouped["reaction_probability"] = grouped["rate"] * total_sim_time_s
    
    total_rxn_prob = grouped["reaction_probability"].sum()
    grouped["selection_probability"] = (
        grouped["reaction_probability"] / total_rxn_prob if total_rxn_prob > 0 else 0.0
    )
    
    grouped["file_id"] = idx
    grouped["system"] = fp.relative_to(run_dir).parts[0]
    grouped["reaction_type"] = grouped["id_diff"].apply(lambda d: f"1-{int(d) + 1}")
    full_dfs.append(grouped)

df_sim_full = pd.concat(full_dfs, ignore_index=True)
df_sim_full["origin"] = "Simulation"

 #%%
# Load Experimental Data
df_exp = get_experimental_dataframe(TARGET_TEMP_K)

# Merge simulation and experimental datasets for comparisons
df_compare = pd.concat(
    [
        df_sim_full[["reaction_type", "rate", "reaction_probability", "origin"]],
        df_exp[["reaction_type", "rate", "reaction_probability", "origin"]],
    ],
    ignore_index=True,
)

reaction_order = sorted(df_compare["reaction_type"].unique())
#%%
# --------------------------------------------------------------------------
# Plot 1: Absolute Rate (Log Scale) - Simulation vs Experiment
# --------------------------------------------------------------------------
plt.figure(figsize=(9, 6))
sns.boxplot(
    data=df_compare,
    x="reaction_type",
    y="rate",
    hue="origin",
    order=reaction_order,
    palette={"Simulation": "#2b5c8f", "Experiment": "#d95f02"},
    width=1,  
    boxprops={"alpha": 0.8}, 
    showfliers=False,
)

# Overlay individual data points
sns.stripplot(
    data=df_compare,
    x="reaction_type",
    y="rate",
    hue="origin",
    order=reaction_order,
    dodge=True,
    color="black",
    jitter=0.15,
    size=5,
    alpha=0.7,
    legend=False,
)

# Annotate "n.a" at y=0.1 for missing reaction types
na_reactions = ["1-2", "1-7", "1-8"]
for rxn in na_reactions:
    if rxn in reaction_order:
        x_pos = reaction_order.index(rxn)
        plt.text(x_pos, 0.1, "n.a", ha="center", va="center", fontsize=10)

plt.yscale("log")
# plt.ylim(0,500)
plt.title(
    f"Absolute Reaction Rates at {TARGET_TEMP_K:.0f} K: Simulation vs Experiment",
    fontsize=13,
)
plt.xlabel("HAT Reaction Type", fontsize=11)
plt.ylabel("Reaction Rate s⁻¹ (log scale)", fontsize=11)
plt.legend(title="Data Source")
plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(output_dir / f"1_absolute_rates_{run_name}.png", dpi=300)
plt.show()


#%% prepare fig 2
# 1. Simulation: Auswahlwahrscheinlichkeiten pro Trajektorie/Datei berechnen
df_sim_prob_final = df_sim_full[["file_id", "system", "reaction_type", "selection_probability"]].copy()
df_sim_prob_final["origin"] = "Simulation"

# --- Concise Creation of Modified Simulation Data ---
zero_types = ["1-2", "1-3", "1-7", "1-8"]

df_sim_masked = df_sim_prob_final.copy().assign(
    selection_probability=lambda df: df["selection_probability"].where(~df["reaction_type"].isin(zero_types), 0.0),
    origin="Simulation (masked)"
)
df_sim_masked["selection_probability"] /= df_sim_masked.groupby("file_id")["selection_probability"].transform("sum")

# 2. Experiment: Mittelwert der Raten über alle Arbeiten bilden & auf Summe = 1.0 normalisieren
exp_mean_rates = (
    df_exp.groupby("reaction_type", as_index=False)["rate"]
    .mean()
    .rename(columns={"rate": "mean_rate"})
)

tot_exp_rate = exp_mean_rates["mean_rate"].sum()
exp_mean_rates["selection_probability"] = (
    exp_mean_rates["mean_rate"] / tot_exp_rate if tot_exp_rate > 0 else 0.0
)
exp_mean_rates["origin"] = "Experiment"

# 3. Combine all three datasets
df_prob_combined = pd.concat(
    [
        df_sim_prob_final[["reaction_type", "selection_probability", "origin"]],
        df_sim_masked[["reaction_type", "selection_probability", "origin"]],
        exp_mean_rates[["reaction_type", "selection_probability", "origin"]],
    ],
    ignore_index=True,
)

# Define ordering and colors
all_reaction_types = sorted(df_prob_combined["reaction_type"].unique())
origin_order = ["Simulation", "Simulation (masked)", "Experiment"]
palette = {
    "Simulation": "#2b5c8f",
    "Simulation (masked)": "#41b6c4",
    "Experiment": "#d95f02",
}
#%%
# --- Brier score divergence, same definition as the paper (fig_3c_and_S2.py + util.brier_score_divergence):
# per simulation run, SUM over the molecule's shift types of (p_sim - p_exp)^2, then MEAN over runs,
# separately per molecule. Summing (not averaging) over shift types makes the score independent of how
# many zero-probability shift types are included.
exp_means = exp_mean_rates.set_index("reaction_type")["selection_probability"]

N_CARBONS = {"propyl": 3, "butyl": 4, "pentyl": 5, "hexyl": 6, "heptyl": 7, "octyl": 8}

# Paper values (fig. 3c), recomputed with this definition from output/n-alkyl_examples/fig_3c.xlsx.
# The paper scored shifts 1-2..1-7 for both molecules (octyl's 1-8 omitted, p ~ 5e-4).
BRIER_PAPER = {"heptyl": 0.0445, "octyl": 0.0200}


def shift_types(system: str) -> list[str]:
    """All possible 1-k H-shifts of an n-alkyl radical with the radical on the terminal carbon."""
    return [f"1-{k}" for k in range(2, N_CARBONS[system] + 1)]


def calc_brier_score(df_sim: pd.DataFrame, system: str) -> float:
    shifts = shift_types(system)
    p_exp = exp_means.reindex(shifts, fill_value=0.0)
    p_sim = (
        df_sim[df_sim["system"] == system]
        .pivot_table(index="file_id", columns="reaction_type", values="selection_probability", aggfunc="sum")
        .reindex(columns=shifts)
        .fillna(0.0)
    )
    return float(((p_sim - p_exp) ** 2).sum(axis=1).mean())


def calc_brier_uniform(system: str) -> float:
    shifts = shift_types(system)
    p_exp = exp_means.reindex(shifts, fill_value=0.0)
    return float(((1.0 / len(shifts) - p_exp) ** 2).sum())


brier = {
    system: dict(
        raw=calc_brier_score(df_sim_prob_final, system),
        masked=calc_brier_score(df_sim_masked, system),
        paper=BRIER_PAPER.get(system, np.nan),
        uniform=calc_brier_uniform(system),
    )
    for system in systems
}
for system, b in brier.items():
    print(f"Brier {system}: " + ", ".join(f"{k} {v:.4f}" for k, v in b.items()))

#%% Plot figure
fig, ax = plt.subplots(figsize=(11, 6))

# A. Common Barplot
sns.barplot(
    data=df_prob_combined,
    x="reaction_type",
    y="selection_probability",
    hue="origin",
    order=all_reaction_types,
    hue_order=origin_order,
    estimator=np.mean,
    errorbar="se",
    palette=palette,
    alpha=0.85,
)

# B. Individual scatter points for both simulation runs
df_sims_combined = pd.concat([df_sim_prob_final, df_sim_masked], ignore_index=True)

sns.stripplot(
    data=df_sims_combined,
    x="reaction_type",
    y="selection_probability",
    hue="origin",
    order=all_reaction_types,
    hue_order=origin_order,
    dodge=True,
    color="black",
    jitter=0.12,
    size=4.0,
    alpha=0.7,
    legend=False,
)

# Diagramm-Formatierung
plt.title(f"HAT Selection Probability at {TARGET_TEMP_K:.0f} K", fontsize=13)
plt.xlabel("HAT Reaction Type", fontsize=11)
plt.ylabel("Selection Probability (linear scale)", fontsize=11)
plt.ylim(0, 1.05)
plt.grid(axis="y", linestyle="--", alpha=0.5)
# Create main legend
leg = plt.legend(
    title="Data Source",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0.,
)

# Format Brier Score text cleanly with line breaks
brier_text = "Brier score vs. Exp\n(sum over shifts, mean over runs):"
for system, b in brier.items():
    brier_text += (
        f"\n{system.capitalize()} ({len(shift_types(system))} shifts)\n"
        f"• Sim raw: {b['raw']:.4f}\n"
        f"• Sim masked: {b['masked']:.4f}\n"
        f"• Paper (fig. 3c): {b['paper']:.4f}\n"
        f"• Uniform: {b['uniform']:.4f}"
    )

# Add text box using fixed axes coordinates aligned under x=1.02
plt.text(
    1.02, 0.65, brier_text,
    transform=ax.transAxes,
    fontsize=9,
    verticalalignment="top",
    linespacing=1.4,
    bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="0.8", alpha=0.9)
)

# Crucial: Specify bbox_extra_artists so tight_layout accounts for the text box on the right
plt.tight_layout()
plt.savefig(output_dir / f"2_selection_probability_{run_name}.png", dpi=300)

#%% Subsampling execution for Figure 3
# Subsample directly from the pre-parsed DataFrames stored in memory
import math

def round_to_nearest_magnitude(val: int) -> int:
    """Rounds a number to its most significant digit (e.g., 790470 -> 800000, 780 -> 800)."""
    if val <= 0:
        return val
    order = 10 ** (len(str(val)) - 1)
    return int(round(val / order) * order)

processed_dfs = [
    process_subsampling(df_raw, idx, *traj_info[fp])
    for idx, (df_raw, fp) in enumerate(zip(raw_parsed_dfs, file_paths))
]

# Concatenate into single dataframe for Plot 3
df_sim = pd.concat(processed_dfs, ignore_index=True)
df_sim["num_frames_rounded"] = df_sim["num_frames"].apply(round_to_nearest_magnitude)
#%%
# --------------------------------------------------------------------------
# Plot 3: Downsampling Effect (Rates vs Frames using 'cividis' palette)
# --------------------------------------------------------------------------
plt.figure(figsize=(12, 6))

frame_order = sorted(df_sim["num_frames_rounded"].unique())
diff_order = sorted(df_sim["reaction_type"].unique())

sns.barplot(
    data=df_sim,
    x="num_frames_rounded",
    y="rate",
    hue="reaction_type",
    order=frame_order,
    hue_order=diff_order,
    estimator=np.mean,
    errorbar="se",
    palette="cividis",
    alpha=0.85,
)

sns.stripplot(
    data=df_sim,
    x="num_frames_rounded",
    y="rate",
    hue="reaction_type",
    order=frame_order,
    hue_order=diff_order,
    dodge=True,
    color="black",
    jitter=0.1,
    size=3.5,
    alpha=0.8,
    legend=False,
)

plt.yscale("log")
plt.title("Effect of Downsampling Resolution on HAT Reaction Rates", fontsize=13)
plt.xlabel("Number of Downsampled Frames", fontsize=11)
plt.ylabel("Mean Reaction Rate s⁻¹ (log scale)", fontsize=11)
plt.legend(title="HAT Reaction", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(output_dir / f"3_rates_vs_frames_cividis_{run_name}.png", dpi=300)
plt.show()

print("Analysis complete. Generated plots successfully saved.")
# %%
#convert s to fs in dfs

for curr_df in raw_parsed_dfs:
    #curr_df['timespans'] = 1e15 * curr_df['timespans']
    #curr_df['timespans'] = curr_df['timespans'].astype(int)
    #curr_df['timespans'] = round(curr_df['timespans'],ndigits=3)
    print(curr_df['timespans'].max())
    print(curr_df['timespans'].min())
# %%
