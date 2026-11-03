"""
03_diet_constraints.py

Apply diets as exchange-reaction bounds to a genome-scale metabolic model and
compare what each diet allows the model to achieve.

For each diet file:
    1. close all uptakes (exchange lower bounds set to 0),
    2. open the diet's uptakes at the listed rates, plus always-open
       inorganic exchanges (water, protons, phosphate, oxygen),
    3. maximise the objective (default: ATP hydrolysis, DM_atp_c_ in Recon3D),
    4. report the reduced cost of each diet uptake, i.e. how much the
       objective would improve per extra unit of that nutrient. Nonzero
       reduced costs identify the nutrients that limit the objective.

Usage:
    # human model (download Recon3D.xml from BiGG first, see README)
    python 03_diet_constraints.py --model Recon3D.xml \
        --diets diets/baseline_diet.csv diets/dairy_enriched_diet.csv

    # quick test on the bundled E. coli core model
    python 03_diet_constraints.py --model textbook --objective ATPM \
        --diets diets/ecoli_test_diet.csv

Diet CSV columns: exchange_id, uptake, note
    uptake is a positive uptake rate (mmol/gDW/h); it becomes lower_bound = -uptake.

Outputs (results/):
    diet_comparison.csv   objective value per diet
    diet_reduced_costs.csv reduced cost of every diet uptake, per diet
    diet_comparison.png
"""

import argparse
from pathlib import Path

import cobra
import matplotlib.pyplot as plt
import pandas as pd

OUT = Path("results")
OUT.mkdir(exist_ok=True)

ALWAYS_OPEN = ["EX_h2o_e", "EX_h_e", "EX_pi_e"]


def load_model(path):
    if path == "textbook":
        return cobra.io.load_model("textbook")
    p = Path(path)
    if p.suffix == ".json":
        return cobra.io.load_json_model(p)
    if p.suffix == ".mat":
        return cobra.io.load_matlab_model(p)
    return cobra.io.read_sbml_model(p)


def apply_diet(model, diet, o2_uptake):
    """Close all uptakes, then open the diet. Returns the diet IDs applied."""
    for ex in model.exchanges:
        ex.lower_bound = 0.0
    for rid in ALWAYS_OPEN:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = -1000.0
    if "EX_o2_e" in model.reactions:
        model.reactions.EX_o2_e.lower_bound = -o2_uptake

    applied, missing = [], []
    for row in diet.itertuples():
        if row.exchange_id in model.reactions:
            model.reactions.get_by_id(row.exchange_id).lower_bound = -float(row.uptake)
            applied.append(row.exchange_id)
        else:
            missing.append(row.exchange_id)
    if missing:
        print(f"  warning: not in model, skipped: {', '.join(missing)}")
    return applied


def run_diet(model, diet_file, objective, o2_uptake):
    diet = pd.read_csv(diet_file, comment="#")
    name = Path(diet_file).stem
    with model:
        model.objective = objective
        model.objective_direction = "max"
        applied = apply_diet(model, diet, o2_uptake)
        sol = model.optimize()
        value = sol.objective_value if sol.status == "optimal" else float("nan")
        rc = (sol.reduced_costs[applied] if sol.status == "optimal"
              else pd.Series(index=applied, dtype=float))
    print(f"  {name}: status={sol.status}, objective={value:.4f}")
    rc_df = pd.DataFrame({"diet": name, "exchange_id": applied,
                          "reduced_cost": rc.values})
    return {"diet": name, "objective": objective, "value": value,
            "status": sol.status}, rc_df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True,
                    help="SBML/JSON/MAT model file, or 'textbook'")
    ap.add_argument("--diets", nargs="+", required=True, help="diet CSV files")
    ap.add_argument("--objective", default="DM_atp_c_",
                    help="reaction to maximise (default: DM_atp_c_)")
    ap.add_argument("--o2", type=float, default=1000.0,
                    help="maximum O2 uptake (default: 1000, unconstrained)")
    args = ap.parse_args()

    print(f"Loading {args.model} ...")
    model = load_model(args.model)
    print(f"  {len(model.reactions)} reactions, {len(model.metabolites)} "
          f"metabolites, {len(model.exchanges)} exchanges")
    if args.objective not in model.reactions:
        hits = [r.id for r in model.reactions if "atp" in r.id.lower()][:15]
        raise SystemExit(f"Objective {args.objective} not found. "
                         f"ATP-related reactions include: {hits}")

    summary, rcs = [], []
    for f in args.diets:
        s, rc = run_diet(model, f, args.objective, args.o2)
        summary.append(s)
        rcs.append(rc)

    summary = pd.DataFrame(summary)
    rcs = pd.concat(rcs, ignore_index=True)
    summary.to_csv(OUT / "diet_comparison.csv", index=False)
    rcs.to_csv(OUT / "diet_reduced_costs.csv", index=False)

    print("\nObjective by diet:")
    print(summary.round(4).to_string(index=False))
    limiting = rcs[rcs["reduced_cost"].abs() > 1e-9]
    print("\nLimiting nutrients (nonzero reduced cost):")
    print(limiting.round(4).to_string(index=False) if len(limiting)
          else "  none")

    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    ax.bar(summary["diet"], summary["value"], color="C0",
           label=f"max {args.objective}")
    ax.set_ylabel("Objective flux (mmol/gDW/h)")
    ax.legend(fontsize=8)
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(OUT / "diet_comparison.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
