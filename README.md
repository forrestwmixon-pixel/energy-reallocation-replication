# Electricity supply reallocation during wind and solar expansion

Replication materials for **System Scale and Trade in Electricity Supply Reallocation during Wind and Solar Expansion: An accounting panel analysis of 199 electricity systems from 2000 to 2024**.

Author: Biao Mu, Nanjing University of Posts and Telecommunications.
Manuscript: Energy Reports, EGYR-D-26-03711.
Replication release: 1.1.0, 2026-11-01.

Public repository: https://github.com/forrestwmixon-pixel/energy-reallocation-replication

Download `EGYR-D-26-03711_Public_Replication_v1.1.0.zip` from this repository and extract it before running the commands below.

## What is included

- `data/owid-energy-data.xlsx`: the unchanged retained OWID source workbook, including its metadata sheet.
- `data/processed_panel.csv`: the full derived panel.
- `code/reanalyse.py`: source extraction, all regression specifications, 9,999-draw wild cluster bootstrap, result tables, and figures.
- `code/practical_checks.py`: realized demand regimes, initial scale restrictions, partial exposure concentration, normalized influence exclusions and within-equation technology contrasts.
- `code/verify_results.py`: an independent full-design OLS/covariance check, a directly re-estimated bootstrap draw, and common-sample Frisch–Waugh–Lovell and adding-up checks.
- `results/`: every full-precision model and sensitivity table, sample flow, system list, source metadata, and verification outputs.
- `figures/`: 300-dpi PNG and vector PDF versions of the four main figures and one supplementary figure.
- `SHA256SUMS.txt`: SHA-256 checksums for the release contents.

This public package excludes correspondence, reviewer reports, and manuscript drafts. It contains the numerical replication material only.

## Reproduce

Use Python 3.12 and the package versions in `requirements.txt`. From the package root:

```bash
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/reanalyse.py --source data/owid-energy-data.xlsx
OPENBLAS_NUM_THREADS=1 python code/verify_results.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/practical_checks.py
```

No live data download or API account is needed. Output files are written into `data`, `results`, and `figures`. The full workflow was rerun from the unchanged workbook in a separate directory before this release. All 31 numerical CSV outputs were byte-identical to the revision outputs. The direct full-design verification also passed. The evidence is in `results/fresh_reproduction_check.json` and `results/verification.json`.

The scripts were reconstructed and verified during revision; they are not represented as recovered original historical scripts.

## Source snapshot

Source: Our World in Data Energy Dataset, integrating the indicator sources identified in the workbook metadata, including Ember and the Energy Institute.

- OWID repository and documentation: https://github.com/owid/energy-data
- Ember source and methodology entry: https://ember-energy.org/data/yearly-electricity-data/
- Retained XLSX SHA-256: `370d84b4c3d09bc2cd3e2996a967b85883625b8e99a982386507281d726559ce`

Later live downloads can differ. The retained source is supplied to avoid that version drift. Missing entries are not filled by the analysis. Source estimates and precision are retained rather than independently revalidated. The old submission checksum could not be verified against a retained source and is superseded by the checksum above.

## Sample and methods

The principal sample contains 4,726 consecutive-year differences from 199 electricity systems over 2001–2024. Coverage eligibility uses wind, solar, fossil generation and demand, with at least 20 joint years and some positive recorded generation. Carbon intensity is not a regression eligibility variable. The component allocation sample contains 4,656 observations and 196 systems; normalization additionally requires positive lagged generation.

Conditional and unconditional regressions use common samples where the coefficients are compared. All models include system and year fixed effects. Conventional inference uses system-clustered CR1 covariance counting the full-rank design and a t(G−1) reference. Wild inference uses 9,999 unrestricted Rademacher cluster bootstrap-t draws, re-estimating all fixed effects, with seed 20261101. Every normalized change, including demand, uses the same positive lagged-generation denominator. The small retained-file discrepancy is kept explicitly instead of adjusting source components to force closure.

## Main numerical output files

- `T01_main.csv`: principal conditional and unconditional models.
- `T02_scale.csv`: scale groups, influential-system exclusions and normalized principal model.
- `T03_allocation.csv`: complete common-sample conditional and unconditional allocation, including all demand coefficients.
- `T04_normalized.csv`: complete normalized allocations overall and by scale group.
- `T05_timing.csv`: isolated and joint lead/lag models, with and without demand, and cumulative tests.
- `T06_technology.csv`: wind and solar jointly, with coal, gas, oil/other fossil, hydro, nuclear, other renewables and net imports.
- `S05_wild_bootstrap.csv`: bootstrap intervals.
- `S09_leave_one_out.csv`: all 199 system omissions.
- `S10_zero_exposure.csv`, `S11_sensitivity.csv`, `S12_composition_groups.csv`: exposure support and additional sensitivity/composition analyses.

The remaining S-prefixed files contain source roles, sample flow, coverage, balance residuals, demand decomposition, normalized unconditional results, persistence, descriptive trajectories and endpoint checks.

## Additional checks in version 1.1.0

Normalization changes which systems carry exposure variation. Falkland Islands account for 42.35% of residualized normalized wind–solar exposure sum of squares in the common sample. Omitting that system leaves the fossil and net-import coefficients at −0.362 and −0.641; omitting the five highest exposure-concentration systems gives −0.448 and −0.541. Thus the fitted trade margin persists, while normalization does not create an equally weighted typical-system estimate.

S18 reports complete normalized component allocations by realized demand growth versus flat/declining demand and by 2000–2004 mean generation floors of 0, 0.1, 1 and 5 TWh (at least three baseline years). S19 reports partial exposure variation by system, after controlling for demand and fixed effects. S20 tests wind minus solar within each joint outcome equation, using the full clustered covariance. S21 reports full normalized influence exclusions. Four S22 files report 9,999-draw wild bootstrap intervals for fossil generation and imports under the two normalized exclusions. practical_verification.json records independent direct-design checks of the additional specifications. All additional checks are exploratory and descriptive; realized demand regimes are not exogenous treatment groups.

## Attribution and reuse

Cite Biao Mu and this replication release when using the analysis. `CITATION.cff` provides the author and release metadata. The software licence in `LICENSE-CODE.txt` applies only to the analysis scripts. Underlying and processed source-derived data retain the original data providers’ applicable terms; this release does not relicense third-party data. The full original indicator attribution is preserved in the source metadata and `results/S03_source_metadata.csv`.
