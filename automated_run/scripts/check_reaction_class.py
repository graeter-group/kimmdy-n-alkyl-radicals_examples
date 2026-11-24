"""Check the HAT reaction class stored in KIMMDY cut-out npz files against the topology.

KIMMDY-HAT (generic_sp3 capping) writes one ``<hash>.npz`` per candidate HAT into
``sim_XX/HAT_000/1_hat_reaction/se/``. Its meta dict contains

    indices         (H id, acceptor radical id, donor heavy-atom id), 1-based global ids
    reaction_class  "X1-n", from a shortest-path search over the *capped cut-out* bond graph
                    ("intermolecular" if donor and acceptor are not connected in the cut-out)

This script recomputes the class from the full-molecule bond graph in the system's
``topol.top`` ([ bonds ] section) and reports every disagreement. It also checks that
the transferred H really is bonded to the stated donor, and compares with the
``1-{|id_rad - id_donor| + 1}`` convention (valid for C1..Cn numbering of the chain).

The se/ directory only survives a successful run if kimmdy.yml has
``reactions: hat_reaction: keep_structures: true``. With the default (false) the
files live in a temporary directory and are deleted.

Usage (from automated_run/):
    python scripts/check_reaction_class.py                 # newest run in runs/
    python scripts/check_reaction_class.py runs/<RUN_NAME> --systems octyl --max-per-sim 2000
    python scripts/check_reaction_class.py --csv mismatches.csv
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:  # tqdm is optional
    def tqdm(it, **_kwargs):
        return it


AUTOMATED_RUN_DIR = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------- topology
def parse_topology(top_path: Path) -> tuple[dict[int, str], dict[int, set[int]]]:
    """Return atom id -> atom name and the bond adjacency from a flat GROMACS .top.

    Only reads [ atoms ] and [ bonds ] of the topology file itself (the n-alkyl
    systems hold a single moleculetype without #include'd molecule itps).
    """
    names: dict[int, str] = {}
    graph: dict[int, set[int]] = defaultdict(set)
    section = None
    for raw in top_path.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        m = re.match(r"\[\s*(\w+)\s*\]", line)
        if m:
            section = m.group(1).lower()
            continue
        fields = line.split()
        if section == "atoms":
            names[int(fields[0])] = fields[4]
        elif section == "bonds":
            a, b = int(fields[0]), int(fields[1])
            graph[a].add(b)
            graph[b].add(a)
    if not names or not graph:
        raise ValueError(f"No [ atoms ] or [ bonds ] found in {top_path}")
    return names, dict(graph)


def bond_distance(graph: dict[int, set[int]], start: int, end: int) -> int | None:
    """Number of bonds on the shortest path start -> end, None if not connected."""
    if start == end:
        return 0
    seen = {start}
    queue = deque([(start, 0)])
    while queue:
        node, dist = queue.popleft()
        for nb in graph.get(node, ()):
            if nb == end:
                return dist + 1
            if nb not in seen:
                seen.add(nb)
                queue.append((nb, dist + 1))
    return None


def topology_class(graph: dict[int, set[int]], donor: int, acceptor: int) -> str:
    """HAT class as KIMMDY defines it: n bonds between donor and acceptor -> '1-{n+1}'."""
    dist = bond_distance(graph, donor, acceptor)
    return "intermolecular" if dist is None else f"1-{dist + 1}"


def normalize_class(cls) -> str:
    """'X1-5' -> '1-5'; anything else unchanged."""
    cls = str(cls)
    return cls[1:] if re.fullmatch(r"X1-\d+", cls) else cls


# --------------------------------------------------------------------------- npz
def load_meta(path: Path) -> dict:
    """Load the meta dict written by np.savez(path, meta) (same as kimmdy_hat)."""
    with np.load(path, allow_pickle=True) as npz:
        data = npz[npz.files[0]]
    if isinstance(data, dict):
        return data
    data = np.asarray(data, dtype=object)
    item = data.item() if data.shape == () else data.reshape(-1)[0]
    return item if isinstance(item, dict) else item.item()


def iter_npz(sim_dir: Path, max_files: int | None):
    """Yield npz paths from sim_XX/HAT_*/*hat_reaction/se/ without globbing whole trees."""
    n = 0
    for hat_dir in sorted(sim_dir.glob("HAT_*")):
        if not hat_dir.is_dir():
            continue
        for task_dir in sorted(hat_dir.glob("*hat_reaction")):
            se_dir = task_dir / "se"
            if not se_dir.is_dir():
                continue
            with os.scandir(se_dir) as it:
                for entry in it:
                    if entry.name.endswith(".npz"):
                        yield Path(entry.path)
                        n += 1
                        if max_files is not None and n >= max_files:
                            return


# --------------------------------------------------------------------------- main
def newest_run(runs_dir: Path) -> Path:
    runs = [p for p in runs_dir.iterdir() if p.is_dir()]
    if not runs:
        raise FileNotFoundError(f"No run directories in {runs_dir}")
    # run names end in _YYYYMMDD_HHMMSS; fall back to mtime
    def key(p: Path):
        m = re.search(r"(\d{8}_\d{6})$", p.name)
        return (m.group(1) if m else "", p.stat().st_mtime)
    return max(runs, key=key)


def check_system(sys_dir: Path, max_per_sim: int | None) -> pd.DataFrame:
    top = sys_dir / "topol.top"
    names, graph = parse_topology(top)
    rows = []
    sims = sorted(p for p in sys_dir.glob("sim_*") if p.is_dir())
    for sim in sims:
        for path in tqdm(list(iter_npz(sim, max_per_sim)), desc=f"{sys_dir.name}/{sim.name}", leave=False):
            row = {"system": sys_dir.name, "sim": sim.name, "npz": path.name}
            try:
                meta = load_meta(path)
            except Exception as e:  # corrupt / partial file
                row["error"] = f"load failed: {e}"
                rows.append(row)
                continue
            idx = tuple(int(i) for i in meta.get("indices", ()))
            if len(idx) < 3:
                row["error"] = f"indices has {len(idx)} entries (capping_method={meta.get('capping_method')})"
                rows.append(row)
                continue
            h, rad, donor = idx[:3]
            kimmdy_cls = meta.get("reaction_class", "missing")
            h_partners = sorted(graph.get(h, ()))
            row.update(
                frame=meta.get("frame"),
                h=h,
                acceptor=rad,
                donor=donor,
                h_name=names.get(h),
                acceptor_name=names.get(rad),
                donor_name=names.get(donor),
                h_bonded_to=h_partners,
                kimmdy_class=str(kimmdy_cls),
                kimmdy_class_norm=normalize_class(kimmdy_cls),
                topology_class=topology_class(graph, donor, rad),
                id_diff_class=f"1-{abs(rad - donor) + 1}",
                translation=meta.get("translation"),
            )
            row["h_on_donor"] = donor in h_partners
            row["class_match"] = row["kimmdy_class_norm"] == row["topology_class"]
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", nargs="?", type=Path,
                        help="run directory (default: newest in automated_run/runs/)")
    parser.add_argument("--systems", nargs="*", help="system subdirs to check (default: all with topol.top)")
    parser.add_argument("--max-per-sim", type=int, default=None,
                        help="read at most this many npz files per sim_XX (limits sshfs traffic)")
    parser.add_argument("--csv", type=Path, help="write all checked rows to this CSV")
    parser.add_argument("--show", type=int, default=20, help="number of mismatch rows to print")
    args = parser.parse_args()

    run_dir = args.run_dir or newest_run(AUTOMATED_RUN_DIR / "runs")
    print(f"Run: {run_dir}")
    systems = args.systems or sorted(p.name for p in run_dir.iterdir() if (p / "topol.top").is_file())

    frames = []
    for system in systems:
        sys_dir = run_dir / system
        if not (sys_dir / "topol.top").is_file():
            print(f"  {system}: no topol.top, skipped")
            continue
        df = check_system(sys_dir, args.max_per_sim)
        print(f"  {system}: {len(df)} npz files")
        frames.append(df)

    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if df.empty:
        print("\nNo npz files found. The se/ directories are only kept with "
              "'keep_structures: true' under reactions.hat_reaction in kimmdy.yml.")
        return 1

    if "error" in df:
        errors = df[df["error"].notna()]
        if len(errors):
            print(f"\n{len(errors)} files could not be checked:")
            print(errors[["system", "sim", "npz", "error"]].head(args.show).to_string(index=False))
        df = df[df["error"].isna()]

    pd.set_option("display.width", 200)
    print("\nKIMMDY reaction_class (rows) vs. topology class (columns):")
    for system, sub in df.groupby("system"):
        print(f"\n[{system}]")
        print(pd.crosstab(sub["kimmdy_class"], sub["topology_class"], margins=True).to_string())

    n = len(df)
    n_cls = int((~df["class_match"]).sum())
    n_h = int((~df["h_on_donor"]).sum())
    n_iddiff = int((df["id_diff_class"] != df["topology_class"]).sum())
    print(f"\nChecked {n} cut-outs")
    print(f"  reaction_class != topology class : {n_cls}")
    print(f"  H not bonded to stated donor     : {n_h}")
    print(f"  id-diff convention != topology   : {n_iddiff}")

    bad = df[~df["class_match"] | ~df["h_on_donor"]]
    if len(bad):
        cols = ["system", "sim", "npz", "frame", "h", "acceptor", "donor", "h_bonded_to",
                "kimmdy_class", "topology_class", "h_on_donor"]
        print(f"\nFirst {min(args.show, len(bad))} problem rows:")
        print(bad[cols].head(args.show).to_string(index=False))
        print("\nMismatch counts (kimmdy -> topology):")
        print(Counter(zip(bad["kimmdy_class"], bad["topology_class"])).most_common())

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"\nWrote {args.csv}")
    return 0 if (n_cls == 0 and n_h == 0) else 2


if __name__ == "__main__":
    sys.exit(main())
