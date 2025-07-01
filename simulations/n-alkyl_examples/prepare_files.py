import numpy as np
from pathlib import Path
import shutil

cwd = Path(__file__).parent

T = 500

n_unique = 100

total_time = 0.75

systems_root_dir = cwd / "systems" / f"{T}K_n{n_unique}"

seeds = np.arange(20)

systems = {
    "propyl": {"radicals: 'X'": "radicals: '3'"},
    "butyl": {"radicals: 'X'": "radicals: '4'"},
    "pentyl": {"radicals: 'X'": "radicals: '5'"},
    "hexyl": {"radicals: 'X'": "radicals: '6'"},
    "heptyl": {"radicals: 'X'": "radicals: '7'"},
    "octyl": {"radicals: 'X'": "radicals: '8'"},
}

base_files = {
    "BASE__kimmdy.yml",
    "BASE__sd.mdp",
    "BASE__sd_slowgrowth.mdp",
    "BASE__submit.sh",
}
base_files = {cwd / "templates" / x for x in base_files}


def replace_in_file(filename, line_patterns):
    with open(filename, 'r') as file:
        lines = file.readlines()

    with open(filename, 'w') as file:
        for line in lines:
            line_stripped = line.strip()

            if line_stripped in line_patterns:
                print(f"Replace {line_stripped} with {line_patterns[line_stripped]}")
                file.write(f"{line_patterns[line_stripped]}\n")
            else:
                file.write(line)


for system in systems.keys():
    print(f"System: {system}")
    wrk_dir = systems_root_dir / system

    slurm_dir = wrk_dir / "slurm"
    slurm_dir.mkdir(exist_ok=True, parents=True)

    # Copy topology files
    for f in (cwd / "templates" / "assets" / f"{system}").glob("*"):
        shutil.copy(f, wrk_dir / f.name)

    for seed in seeds:
        for base_file in base_files:
            if base_file.name != "BASE__submit.sh":
                run_file = wrk_dir / f'{base_file.stem.split("BASE__")[1]}__s{seed}{base_file.suffix}'
            else:
                run_file = wrk_dir / f'{base_file.stem.split("BASE__")[1]}{base_file.suffix}'

            print(f"Copy {base_file} to {run_file}")
            shutil.copy(base_file, run_file)

            line_patterns = {}
            if base_file.name == "BASE__kimmdy.yml":
                line_patterns.update(systems[system])
                line_patterns.update({
                    "name: 'HAT_000'": f"name: '{system}_hat_s{seed:03d}__000'",
                    "top: 'X'": f"top: '{system}R.top'",
                    "gro: 'X'": f"gro: '{system}R_box.gro'",
                    "mdp: 'sd.mdp'": f"    mdp: 'sd__s{seed}.mdp'",
                    "mdp: 'sd_slowgrowth.mdp'": f"    mdp: 'sd_slowgrowth__s{seed}.mdp'",
                    "temperature: X": f"      temperature: {T}",
                    "n_unique: 500": f"    n_unique: {n_unique}"
                })
                replace_in_file(run_file, line_patterns)
            elif base_file.name == "BASE__sd.mdp" or base_file.name == "BASE__sd_slowgrowth.mdp":
                line_patterns.update({
                    "gen_seed = X": f"gen_seed = {seed}",
                    "ld-seed = X":  f"ld-seed = {seed}",
                    "ref-t                    = X": f"ref-t                    = {T}",
                })
                replace_in_file(run_file, line_patterns)
            elif base_file.name == "BASE__submit.sh":
                line_patterns.update({
                    '#SBATCH --job-name="BASENAME"': f'#SBATCH --job-name="KMDY_{system.upper()}"',
                    '#SBATCH --array=X-Y': f'#SBATCH --array={seeds.min()}-{seeds.max()}%5'
                })
                replace_in_file(run_file, line_patterns)
            else:
                raise ValueError(f"Unrecognized base_file: {base_file}")
