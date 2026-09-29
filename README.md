# Tau-5* ASTEROIDS-style selection

This is an independent reimplementation of fixed-size equal-weight ensemble selection, not the original ASTEROIDS source code. The historical continuously weighted fits remain in [tau5-sec-saxs-analysis](https://github.com/paulrobustelli/tau5-sec-saxs-analysis).

## Current status

The selector and output pipeline are implemented and tested. Full WT/AA selections are waiting for candidate generation and forward calculations; no final fit or convergence claim is made yet.

## Method

The representation follows [Huang et al., JACS 2014, DOI 10.1021/ja502030n](https://doi.org/10.1021/ja502030n): each gene is one conformer, ensemble size stays fixed, all selected conformers have population 1/N, and selection uses mutation and crossover. Here N=500 is the requested size, not a size established by held-out cross-validation.

Our explicit implementation choices are 100 chromosomes, 1,000 generations, and three random seeds per independent candidate pool. Each generation proposes internal mutations from the current population, external mutations from the entire candidate pool, crossover from pairs of parents, and fresh random ensembles. Random four-way tournaments retain offspring; elitism preserves the best solution. All chromosomes contain 500 distinct candidate IDs. Mutation sizes mix 1, 2, 5, and N/20 replacements. These hyperparameters and details are our choices, not a claim to reproduce unavailable authors' code exactly. Runs stop at the generation budget, not a predictor-error tolerance. Further search is required if objective traces or independent runs remain unstable.

The objective is the sum of squared scaled chemical-shift residuals plus the sum of squared SAXS residuals normalized by experimental errors. Carbon scales are 0.5 ppm, N 2.45 ppm, HN 0.49 ppm, and HA 0.25 ppm where measured. These specify the optimization tradeoff; 0.5 ppm is not a claimed predictor accuracy or an acceptance cutoff. Per-nucleus RMSDs and SAXS mean chi-square are reported separately. Fit a single nonnegative SAXS intensity scale per ensemble, with no additive background or chemical-shift offset. There is no entropy penalty and no optimization of individual conformer weights.

## Inputs

GP + AR330–447 (120 residues), WT and W397A/W433A. Model residue = AR residue − 327. Every PDB sequence is checked against its construct. Retain every measured interior shift independently for each construct; missing forward predictions cause an explicit failure. Exclude GP and terminal P447 shifts. Subtract the same SPARTA+ random-coil baseline from experiment and calculation for secondary-shift displays. Report DSSP H and H/G/I helicity, without Ramachandran-basin panels.

Independent candidate groups combine pool_1 + pool_3 or pool_2 + pool_4, each containing 2,000 base structures, with matching helix-length supplements (110 WT or 90 AA per group). Supplements cover lengths 4/6/8/10/12 around sites with D2D helicity ≥10%. IDPConformerGenerator supplies candidates; SPARTA+ and CRYSOL supply predictions. Enrichment frequencies are not physical populations.

## Running

From the original workspace root with the existing prediction cache:

```sh
OPENBLAS_NUM_THREADS=1 work/ensemble_env/bin/python work/asteroids/test_selection.py
OPENBLAS_NUM_THREADS=1 work/ensemble_env/bin/python work/asteroids/run_selection.py WT 1
```

`selection.py` needs NumPy. `run_selection.py` expects experimental inputs under `outputs/ASTEROIDS_setup/` and per-conformer result.json caches under `outputs/Tau5_joint_refinement/{WT,AA}/backcalc/`. Those large caches and licensed predictor programs are not bundled here yet. `wait_and_select.py` waits for exact completed candidate counts before dispatching all four independent groups.

Outputs go to `outputs/Tau5_ASTEROIDS/`: candidate provenance, input fingerprint, equal-weight selected membership, shift and SAXS residual tables, DSSP populations, search traces, and summaries. Best-so-far checkpoints are written every 25 generations. Completed matching runs are skipped; an interrupted seed restarts deterministically from generation zero. A checkpoint is not full optimizer-state continuation.

Validation includes a synthetic target whose optimum is exhaustively enumerated, unique subset cardinality, elitism, the full-pool edge case, and a 500-member end-to-end output smoke test. Real-data fit quality and convergence remain to be assessed.
