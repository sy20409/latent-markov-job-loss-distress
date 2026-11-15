# A Latent Markov Analysis of Job Loss, Continued Unemployment and Re-employment

Seminar paper, *Multivariate Statistics*, University of Trier (2026).
Author: Shinyoung Lee

Examines how job loss, continued unemployment and re-employment change the
probabilities of moving between latent states of psychological distress
(Jahoda's latent deprivation model). The central question is whether unemployment
harms mental health mainly by **triggering deterioration** or by **preventing
recovery**. Distress is treated as a latent state, and a latent Markov model
estimates the transitions between states while correcting for classification error.

## Summary
- **Data**: Understanding Society (UKHLS), Waves 6-10, 9,638 adults aged 25-60
  employed or unemployed at every wave (46,585 person-waves)
- **Measure**: 12 binary GHQ-12 items
- **Methods**: latent class analysis (`poLCA`) to choose three distress states, then a
  latent Markov model with covariate-dependent transitions (`LMest`)
- **Key finding**: job loss raised the odds of moving from low to high distress
  (OR = 4.06); continued unemployment lowered the odds of recovering from high
  distress (OR = 0.15)

## Hypotheses
Derived from Jahoda (1982):
- **H1**: Job loss increases the probability of moving from low to higher distress.
- **H2**: Continued unemployment decreases the probability of recovering from high distress.
- **H3**: Re-employment increases the probability of recovering from moderate or high
  distress compared with remaining unemployed.

## Sample and measures
- Respondents aged 25-60 who were employed or unemployed at **all five waves**
  (retirement, education, family care and long-term sickness are excluded)
- Employment transitions relative to "stayed employed": 417 job losses,
  471 cases of continued unemployment, 448 re-employments
- GHQ-12 items coded 0-0-1-1 (1 = distress); controls are age (per decade, centred
  at 40) and sex
- Waves 6-10 were chosen because they largely precede the COVID-19 pandemic

## Methods
1. **Latent class analysis** on pooled person-waves (2-5 classes), compared by BIC,
   relative entropy and the size of the smallest class. Three states were selected:
   low, moderate and high distress.
2. **Latent Markov models** (`LMest`) with the measurement model and the transitions
   estimated jointly; covariates enter the transitions through a multinomial logit:
   - Model A: no covariates
   - Model B: + employment-transition dummies
   - Model C: + age and sex
3. Odds ratios with 95% Wald confidence intervals; Holm adjustment as the primary
   and Benjamini-Hochberg as the secondary correction for multiple testing.

## Main results
| Finding | Estimate |
|---|---|
| Job loss: low to moderate distress | OR = 2.72 |
| Job loss: low to high distress | OR = 4.06 |
| Continued unemployment: high to low distress (recovery) | OR = 0.15 |
| Re-employment vs. continued unemployment: high to low | OR = 8.09 (95% CI 2.73-23.99) |

Unemployment appears to harm mental health both by triggering deterioration and by
preventing recovery, and re-employment restores the chance of recovery from high
distress. Distress was already elevated before job loss, so the estimates may
overstate its effect.

## Sensitivity analyses
Four latent states, unbalanced panel (at least three waves), lead (timing) test,
duration of unemployment, and GHQ-12 without the two employment-related items.
All tables are in the appendix of the paper.

## Files
- [`Transitions_in_Employment_and_Wellbeing_with_R.pdf`](Transitions_in_Employment_and_Wellbeing_with_R.pdf): full paper
- [`ukhls_employment_lta.R`](ukhls_employment_lta.R): data extraction, sample selection,
  latent class analysis, latent Markov models, sensitivity analyses, tables and figures
- `.gitignore`: keeps the licensed data files out of the repository

## Data
The data are **not included** (UKDS End User Licence does not allow redistribution).
To reproduce the analysis:
1. Register with the [UK Data Service](https://ukdataservice.ac.uk) and download
   *Understanding Society* (SN 6614) in Stata format.
2. Unzip it so that the wave files (for example `f_indresp.dta`) are in
   `UKDA-6614-stata/stata/stata14_se/ukhls/`.
3. Set `project` at the top of the script to the folder that contains `UKDA-6614-stata`.

The sample checks in the script (for example 46,585 person-waves) refer to the
21st edition used in the paper and may differ for other releases.

## Limitations
- Employment status is observed only at annual interviews, so short spells between
  interviews are not captured.
- Observational data: health selection into unemployment cannot be ruled out.

## Tools
R (`haven`, `dplyr`, `poLCA`, `LMest`, `ggplot2`)
