# Tau-5* ASTEROIDS-style ensemble selection

Independent implementation in development for androgen-receptor Tau-5* WT and W397A/W433A (AA). This is not the original authors' ASTEROIDS source code. No completed equal-weight fits are reported yet.

## Planned selection

- Select approximately 500 distinct conformers per ensemble, with exactly equal weight 1/N for each selected member.
- Use fixed-size genetic-algorithm subset selection, including mutation and crossover; do not optimize individual continuous weights.
- Fit all available interior backbone chemical shifts for each construct and SAXS jointly. Minimize shift residuals as closely as the candidate pool permits, with predictor error providing context rather than a tolerance that stops fitting. Report per-nucleus RMSDs and SAXS residuals separately.
- Compare independent candidate pools and repeated selection seeds. A small residual alone does not establish ensemble convergence.

## Inputs and modeling conventions

The construct is GP + AR330–447 (120 residues). Model residues 1–2 are GP; model residue = AR residue − 327. AA contains W397A and W433A, with the remaining sequence unchanged.

Candidate generation targets 2,000 base conformers per independent replicate using IDPConformerGenerator, plus separately tracked helix-length supplements. Supplement windows span 4, 6, 8, 10, and 12 residues at D2D helicity ≥10% sites. These enriched pool frequencies are not physical populations.

SPARTA+ supplies chemical-shift predictions and sequence-specific random-coil references; CRYSOL supplies SAXS predictions. Exclude terminal database-average shifts and unmeasured GP shifts. Use all other measured nuclei separately for each protein. Subtract the same SPARTA+ random-coil reference from experimental and calculated shifts when displaying secondary shifts. Compare helicity using DSSP H and H/G/I assignments with D2D.

Licensed predictor executables are not distributed here. Algorithm settings, objective normalization, and reproducible runs will be documented alongside the implementation before results are reported.

## Comparison baseline

The existing continuous-weight analysis remains in [tau5-sec-saxs-analysis](https://github.com/paulrobustelli/tau5-sec-saxs-analysis), including the [executed comparison notebook](https://github.com/paulrobustelli/tau5-sec-saxs-analysis/blob/main/Chemical_shift_SAXS/Tau5_current_results.ipynb). Its weighted fits are historical comparisons, not ASTEROIDS equal-weight selections.
