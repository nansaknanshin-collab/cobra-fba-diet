"""
02_nutrient_scan.py

How nutrient availability constrains metabolism, using the E. coli core model.

    (a) Oxygen scan: growth and acetate secretion as oxygen uptake falls
        (aerobic to anaerobic; overflow metabolism).
    (b) Carbon source comparison: growth on different single carbon sources
        at the same carbon uptake (60 mmol C/gDW/h).

Outputs (results/):
    oxygen_scan.csv, oxygen_scan.png
    carbon_sources.csv, carbon_sources.png
"""

from pathlib import Path

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path("results")
OUT.mkdir(exist_ok=True)

# exchange reaction -> number of carbon atoms in the substrate
CARBON_SOURCES = {
    "EX_glc__D_e": ("glucose", 6),
    "EX_fru_e": ("fructose", 6),
    "EX_pyr_e": ("pyruvate", 3),
    "EX_succ_e": ("succinate", 4),
    "EX_mal__L_e": ("malate", 4),
    "EX_ac_e": ("acetate", 2),
    "EX_etoh_e": ("ethanol", 2),
}
C_UPTAKE = 60.0  # mmol carbon / gDW / h


def oxygen_scan(model):
    rows = []
    for o2 in np.linspace(0, 25, 26):
        with model:
            model.reactions.EX_o2_e.lower_bound = -o2
            sol = model.optimize()
            ok = sol.status == "optimal"
            rows.append({
                "o2_uptake": o2,
                "growth": sol.objective_value if ok else 0.0,
                "acetate": sol.fluxes["EX_ac_e"] if ok else 0.0,
                "ethanol": sol.fluxes["EX_etoh_e"] if ok else 0.0,
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "oxygen_scan.csv", index=False)

    fig, ax1 = plt.subplots(figsize=(5, 3.5))
    ax1.plot(df["o2_uptake"], df["growth"], "o-", label="growth")
    ax1.set_xlabel("O$_2$ uptake (mmol/gDW/h)")
    ax1.set_ylabel("Growth (1/h)")
    ax2 = ax1.twinx()
    ax2.plot(df["o2_uptake"], df["acetate"], "s--", color="C1",
             label="acetate secretion")
    ax2.plot(df["o2_uptake"], df["ethanol"], "^:", color="C2",
             label="ethanol secretion")
    ax2.set_ylabel("Secretion (mmol/gDW/h)")
    lines = ax1.get_legend_handles_labels()
    lines2 = ax2.get_legend_handles_labels()
    ax1.legend(lines[0] + lines2[0], lines[1] + lines2[1],
               loc="center right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "oxygen_scan.png", dpi=200)
    plt.close(fig)
    return df


def carbon_sources(model):
    rows = []
    for ex_id, (name, n_c) in CARBON_SOURCES.items():
        with model:
            for ex in CARBON_SOURCES:  # close all listed carbon uptakes
                model.reactions.get_by_id(ex).lower_bound = 0.0
            model.reactions.get_by_id(ex_id).lower_bound = -C_UPTAKE / n_c
            sol = model.optimize()
            growth = sol.objective_value if sol.status == "optimal" else 0.0
        rows.append({"source": name, "uptake_mmol": C_UPTAKE / n_c,
                     "growth": growth})
    df = pd.DataFrame(rows).sort_values("growth", ascending=False)
    df.to_csv(OUT / "carbon_sources.csv", index=False)

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.bar(df["source"], df["growth"], color="C0")
    ax.set_ylabel("Growth (1/h)")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(OUT / "carbon_sources.png", dpi=200)
    plt.close(fig)
    return df


def main():
    model = cobra.io.load_model("textbook")
    print("Oxygen scan:")
    print(oxygen_scan(model).round(3).to_string(index=False))
    print(f"\nCarbon sources at {C_UPTAKE:.0f} mmol C/gDW/h:")
    print(carbon_sources(model).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
