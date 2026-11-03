# Constraint-based metabolic modelling with COBRApy: FBA, FVA and diet constraints

A compact, reproducible introduction to constraint-based reconstruction and analysis (COBRA) in Python. It moves from flux balance analysis on a core bacterial model to applying diets as exchange constraints on the human genome-scale model Recon3D.

## Background

Constraint-based modelling represents metabolism by its stoichiometric matrix S. At steady state, internal metabolite concentrations do not change, so the flux vector v must satisfy

    S v = 0,   lb <= v <= ub.

Flux balance analysis (FBA) chooses, among all feasible flux vectors, one that maximises a biological objective (growth, ATP production) by linear programming. Flux variability analysis (FVA) then finds the range each flux can take while the objective stays near its optimum. Nutrient availability enters through the bounds on exchange reactions, which is how a diet is imposed on a human model.

This complements my mechanistic modelling of enzyme kinetics (Michaelis–Menten and product-inhibition ODE and PDE models of polymer degradation). Mechanistic models resolve rates and dynamics for a few reactions; constraint-based models trade kinetic detail for genome-scale coverage.

## Contents

| Script | What it does | Model |
| --- | --- | --- |
| `01_fba_basics.py` | FBA, FVA at 90% of optimum, single-gene deletion screen for essential genes | E. coli core (bundled) |
| `02_nutrient_scan.py` | Oxygen scan (aerobic to anaerobic, overflow metabolism); growth on seven carbon sources at equal carbon uptake | E. coli core (bundled) |
| `03_diet_constraints.py` | Applies diets from CSV files as exchange bounds, maximises an objective, and uses reduced costs to identify limiting nutrients | Recon3D (human) or any SBML/JSON/MAT model |
| `diets/` | Baseline and dairy-enriched illustrative diets (BiGG IDs), plus two E. coli test diets | |

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Run

```bash
python 01_fba_basics.py
python 02_nutrient_scan.py

# test the diet pipeline on the bundled model
python 03_diet_constraints.py --model textbook --objective ATPM \
    --diets diets/ecoli_test_diet.csv diets/ecoli_test_diet_mixed.csv

# human model: download Recon3D (SBML) from BiGG first
#   http://bigg.ucsd.edu/models/Recon3D  ->  Recon3D.xml
python 03_diet_constraints.py --model Recon3D.xml \
    --diets diets/baseline_diet.csv diets/dairy_enriched_diet.csv
```

Outputs (CSV tables and PNG figures) are written to `results/`.

## Results on the E. coli core model

- **FBA:** maximum growth on 10 mmol/gDW/h glucose is 0.874 1/h. Only 8 of 95 fluxes are fixed at 90% of optimal growth, so the network has substantial redundancy.
- **Gene essentiality:** 7 genes are essential in silico (growth below 1% of wild type).
- **Oxygen:** as oxygen uptake falls, growth drops roughly linearly while acetate secretion rises (overflow metabolism). Below about 7 mmol/gDW/h ethanol secretion appears as fermentation takes over (`results/oxygen_scan.png`).
- **Carbon sources:** at equal carbon uptake (60 mmol C/gDW/h), glucose and fructose support the same growth. Pyruvate, succinate, acetate and malate support less, because entering central metabolism downstream of glycolysis costs energy for gluconeogenesis. Ethanol supports the most growth: it is the most reduced substrate, so each carbon delivers more electrons to respiration. With oxygen unconstrained the model can exploit that fully, which is a reminder that FBA predictions are only as realistic as their bounds.
- **Diet pipeline:** adding succinate and glutamine to glucose raises maximum ATP hydrolysis from 175 to 238.75 mmol/gDW/h. All three nutrients have nonzero reduced costs, so each is limiting.

## Interpreting reduced costs

A reduced cost is the change in the objective per unit change in that reaction's bound. Diet uptakes are imposed as negative lower bounds, so a limiting nutrient shows a negative reduced cost; its magnitude is the extra objective per extra mmol of that nutrient. A zero reduced cost means the nutrient is in surplus for this objective.

## Limitations

- The human diets are **illustrative**: rates are in mmol/gDW/h and were chosen to show the method, not taken from dietary intake data. A real analysis would convert intake (for example from a VMH diet or a food composition database) into exchange bounds.
- Maximising ATP is a simple energetic objective. Questions about nutrition and health need objectives tied to specific functions (biosynthesis of a metabolite, tissue-specific tasks) and, for the gut, community models of microbes and host.
- Diet files use BiGG identifiers for Recon3D. Human-GEM uses different reaction IDs, so it needs an ID mapping first.

## Next steps

- Use published VMH diets (for example the Average European Diet) with tissue-specific or whole-body models.
- Add microbial community modelling of the gut microbiome (microbe–host metabolite exchange).
- Reproduce the same analyses in the MATLAB COBRA Toolbox to compare results across the two implementations.
