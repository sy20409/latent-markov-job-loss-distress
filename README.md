# A Latent Markov Analysis of Job Loss, Continued Unemployment and Re-employment

Seminar paper, Multivariate Statistics, University of Trier (2026).
Examines how job loss, continued unemployment and re-employment change the
probabilities of moving between latent states of psychological distress
(Jahoda's latent deprivation model).

## Summary
- Data: Understanding Society (UKHLS), Waves 6-10, 9,638 adults aged 25-60
  employed or unemployed at every wave (46,585 person-waves)
- Measure: 12 binary GHQ-12 items
- Methods: latent class analysis (poLCA) to choose three distress states, then
  a latent Markov model with covariate-dependent transitions (LMest)
- Key finding: job loss raised the odds of moving from low to high distress
  (OR = 4.06); continued unemployment lowered the odds of recovering from high
  distress (OR = 0.15)

## Files
- `<<Transitions in Employment and Wellbeing with R>>.pdf`: full paper
- `well-being code`: data extraction, sample selection, LCA, latent Markov
  models, sensitivity analyses, tables and figures

## Data
The data are not included (UKDS End User Licence). Register with the UK Data
Service and download Understanding Society (SN 6614), then set the folder
paths at the top of the script.

## Tools
R (haven, dplyr, poLCA, LMest)
