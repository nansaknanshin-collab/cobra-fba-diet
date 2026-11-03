"""
01_fba_basics.py

Flux balance analysis (FBA), flux variability analysis (FVA) and in silico
gene essentiality on the E. coli core metabolic model (bundled with COBRApy).

Outputs (results/):
    fba_fluxes.csv        optimal flux distribution
    fva_90pct.csv         flux ranges with growth held at >= 90% of optimum
    gene_essentiality.csv growth after each single-gene deletion
"""

from pathlib import Path

import cobra
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion

OUT = Path("results")
OUT.mkdir(exist_ok=True)


def main():
    model = cobra.io.load_model("textbook")  # E. coli core, 95 reactions
    print(f"Model: {model.id} | {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

    # 1. FBA: maximise biomass under default (aerobic, glucose-limited) bounds
    sol = model.optimize()
    print(f"\nFBA status: {sol.status}")
    print(f"Optimal growth rate: {sol.objective_value:.4f} 1/h")
    print("\nModel summary at optimum:")
    print(model.summary())
    sol.fluxes.rename("flux").to_csv(OUT / "fba_fluxes.csv")

    # 2. FVA: which fluxes are fixed and which are flexible near the optimum?
    fva = flux_variability_analysis(model, fraction_of_optimum=0.9)
    fva["range"] = fva["maximum"] - fva["minimum"]
    fva.sort_values("range", ascending=False).to_csv(OUT / "fva_90pct.csv")
    n_fixed = int((fva["range"].abs() < 1e-6).sum())
    print(f"\nFVA (90% of optimum): {n_fixed} of {len(fva)} reactions have a "
          f"fixed flux; the rest can vary without reducing growth below 90%.")

    # 3. Single-gene deletions: which genes are essential for growth?
    wt = sol.objective_value
    dele = single_gene_deletion(model)
    dele["gene"] = dele["ids"].apply(lambda s: next(iter(s)))
    dele["relative_growth"] = dele["growth"].fillna(0) / wt
    dele = dele[["gene", "growth", "relative_growth", "status"]]
    dele.sort_values("relative_growth").to_csv(
        OUT / "gene_essentiality.csv", index=False)
    essential = dele[dele["relative_growth"] < 0.01]
    print(f"\nEssential genes (growth < 1% of wild type): {len(essential)}")
    print(", ".join(sorted(essential["gene"])))


if __name__ == "__main__":
    main()
