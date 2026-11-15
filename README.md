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
- GHQ-12 items coded 0-0-1-1 (1 = distress); controls
