# =============================================================================
# Transitions in Employment and Well-being
# Latent Markov analysis of GHQ-12 distress states, Understanding Society
# (UKHLS) waves 6-10
# =============================================================================

rm(list = ls())

project <- "C:/Users/LSY/Desktop/Project"
ukhls   <- file.path(project, "UKDA-6614-stata/stata/stata14_se/ukhls")
setwd(project)

for (p in c("haven", "dplyr", "poLCA", "LMest", "ggplot2")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(haven)
library(dplyr)
library(poLCA)   # attaches MASS, which masks dplyr::select
library(LMest)
library(ggplot2)
packageVersion("LMest")

waves    <- c("f", "g", "h", "i", "j")
items    <- paste0("scghq", letters[1:12])
emp_vars <- c("job_loss", "stayed_unemp", "re_employed")
labs3    <- c("Low", "Moderate", "High")
idx      <- c("pidp", "t")
resp     <- as.formula(paste(paste(items, collapse = " + "), "~ NULL"))


# ---- 1. Extract waves 6-10 (f-j) ---------------------------------------------

hhvars     <- "country fihhmnnet1_dv gor_dv hhden01_xw hhden91_xw hhden99_xw hhdeng2_xw hhdenub_xw hhdenui_xw hhdenus_xw hhsize hhtype_dv hidp hsivlw ieqmoecd_dv nkids_dv psu strata tenure_dv urban_dv xpgaslw xpleclw"
indallvars <- "age_dv country ethn_dv ff_ivlolw gor_dv hhsize hidp iviolw mastat_dv nchild_dv pidp pno psnen01_lw psnen01_xw psnen91_lw psnen91_xw psnen99_lw psnen99_xw psneng2_lw psneng2_xw psnenub_lw psnenub_xw psnenui_lw psnenui_xw psnenus_lw psnenus_xw psu racel_dv sex_dv strata urban_dv"
indvars    <- "age_dv country ethn_dv ff_emplw ff_ivlolw ffbrfedlw ftexw gor_dv hhsize hhtype_dv hidp ind5mus_lw ind5mus_xw indbd91_lw indbdub_lw indin01_lw indin01_xw indin91_lw indin91_xw indin99_lw indin99_xw inding2_lw inding2_xw indinub_lw indinub_xw indinui_lw indinui_xw indinus_lw indinus_xw indns91_lw indnsub_lw indpxg2_xw indpxub_lw indpxub_xw indpxui_lw indpxui_xw indpxus_lw indpxus_xw iviolw jbhrs jbstat mastat_dv nbornlw nchild_dv paygu_dv fimnlabnet_dv pidp pno psu racel_dv scghq1_dv scghq2_dv sclfsato sex_dv strata tenure_dv ukborn urban_dv"

prefix_vars   <- function(vars, wave) paste0(wave, "_", strsplit(vars, " ")[[1]])
keep_existing <- function(vars, df) unique(vars[vars %in% names(df)])

extract_wave <- function(wave) {
  hidp  <- paste0(wave, "_hidp")
  hh_v  <- prefix_vars(hhvars, wave)
  all_v <- prefix_vars(indallvars, wave)
  ind_v <- prefix_vars(indvars, wave)

  hhresp <- read_dta(file.path(ukhls, paste0(wave, "_hhresp.dta")))
  hhresp <- dplyr::select(hhresp, all_of(keep_existing(c(hidp, hh_v), hhresp)))

  indall <- right_join(hhresp,
                       read_dta(file.path(ukhls, paste0(wave, "_indall.dta"))),
                       by = hidp)
  indall <- indall[!is.na(indall$pidp), ]
  indall <- dplyr::select(indall, all_of(keep_existing(c("pidp", hidp, hh_v, all_v), indall)))

  indresp <- read_dta(file.path(ukhls, paste0(wave, "_indresp.dta")))
  indresp <- indresp[, !names(indresp) %in% names(indall) | names(indresp) == "pidp"]
  d <- merge(indall, indresp, by = "pidp", all = TRUE)
  d <- dplyr::select(d, all_of(keep_existing(c("pidp", hidp, hh_v, ind_v, all_v), d)))

  names(d)   <- sub(paste0("^", wave, "_"), "", names(d))
  d$wavename <- which(letters == wave)
  zap_labels(d)
}

df_final <- bind_rows(lapply(waves, extract_wave))
df_final <- df_final[order(df_final$pidp), ]
table(df_final$wavename)

saveRDS(df_final, file.path(project, "UKHLS_long_hij.rds"))

ghq_raw <- bind_rows(lapply(waves, function(w) {
  d <- read_dta(file.path(ukhls, paste0(w, "_indresp.dta")),
                col_select = c("pidp", paste0(w, "_", items)))
  names(d)   <- sub(paste0("^", w, "_"), "", names(d))
  d$wavename <- which(letters == w)
  zap_labels(d)
}))


# ---- 2. Analytic samples -----------------------------------------------------
# Aged 25-60, employed (jbstat 1, 2) or unemployed (jbstat 3), observed in at
# least `min_waves` waves. Employment transitions are defined only when the
# previous wave is observed; after_gap flags rows whose previous wave is missing.
# GHQ-12 items: 1 = no distress, 2 = distress (0-0-1-1 scoring).

prepare_sample <- function(min_waves) {
  df_clean <- df_final %>%
    filter(age_dv >= 25 & age_dv <= 60) %>%
    mutate(emp_status = case_when(jbstat %in% c(1, 2) ~ "employed",
                                  jbstat == 3         ~ "unemployed",
                                  TRUE                ~ "other")) %>%
    filter(emp_status != "other")

  keep_ids <- df_clean %>%
    group_by(pidp) %>%
    summarise(n_waves = n_distinct(wavename)) %>%
    filter(n_waves >= min_waves) %>%
    pull(pidp)

  df_trans <- df_clean %>%
    filter(pidp %in% keep_ids) %>%
    arrange(pidp, wavename) %>%
    group_by(pidp) %>%
    mutate(emp_lag   = lag(emp_status),
           wave_lag  = lag(wavename),
           consec    = !is.na(wave_lag) & (wavename - wave_lag == 1),
           after_gap = as.integer(wavename > 6 & !consec),
           transition = case_when(
             !consec ~ NA_character_,
             emp_lag == "employed"   & emp_status == "employed"   ~ "stayed_employed",
             emp_lag == "employed"   & emp_status == "unemployed" ~ "job_loss",
             emp_lag == "unemployed" & emp_status == "employed"   ~ "re_employed",
             emp_lag == "unemployed" & emp_status == "unemployed" ~ "stayed_unemployed")) %>%
    mutate(transition_lag  = lag(transition),
           transition_lead = if_else(lead(wavename) - wavename == 1,
                                     lead(transition), NA_character_)) %>%
    ungroup()

  d <- left_join(df_trans, ghq_raw, by = c("pidp", "wavename"))
  d[items] <- lapply(d[items], function(x) {
    x <- as.numeric(x)
    x[x < 1] <- NA
    ifelse(x >= 3, 2, 1)
  })
  d[complete.cases(d[items]), ]
}

make_lm_data <- function(d) {
  d %>%
    arrange(pidp, wavename) %>%
    mutate(across(all_of(items), ~ . - 1),   # LMest expects categories 0/1
           t            = wavename - 5,
           job_loss     = as.integer(transition %in% "job_loss"),
           stayed_unemp = as.integer(transition %in% "stayed_unemployed"),
           stayed_emp   = as.integer(transition %in% "stayed_employed"),
           re_employed  = as.integer(transition %in% "re_employed"),
           female       = as.integer(sex_dv == 2),
           age_c        = (age_dv - 40) / 10,
           # next wave's transition (lead test); lead_missing = no next wave
           job_loss_lead     = as.integer(transition_lead %in% "job_loss"),
           stayed_unemp_lead = as.integer(transition_lead %in% "stayed_unemployed"),
           re_employed_lead  = as.integer(transition_lead %in% "re_employed"),
           lead_missing      = as.integer(is.na(transition_lead)),
           # duration: second vs. third or later consecutive wave of unemployment
           stayed_unemp_2w = as.integer(transition %in% "stayed_unemployed" &
                                        !(transition_lag %in% "stayed_unemployed")),
           stayed_unemp_3w = as.integer(transition %in% "stayed_unemployed" &
                                        transition_lag %in% "stayed_unemployed")) %>%
    as.data.frame()
}

df_lca <- prepare_sample(min_waves = 5)
nrow(df_lca)
mean(rowSums(df_lca[items] == 2) == as.numeric(df_lca$scghq2_dv))   # 1


# ---- 3. Number of latent states: pooled latent class analysis ----------------

f_lca <- as.formula(paste0("cbind(", paste(items, collapse = ","), ") ~ 1"))

entropy_rel <- function(m) {
  p <- m$posterior
  1 - (-sum(p * log(p + 1e-12))) / (nrow(p) * log(ncol(p)))
}

set.seed(123)
fits_lca <- lapply(2:5, function(k)
  poLCA(f_lca, data = as.data.frame(df_lca), nclass = k,
        nrep = 5, maxiter = 2000, verbose = FALSE))

data.frame(K             = 2:5,
           AIC           = sapply(fits_lca, `[[`, "aic"),
           BIC           = sapply(fits_lca, `[[`, "bic"),
           Entropy       = sapply(fits_lca, entropy_rel),
           Min_class_pct = sapply(fits_lca, function(m) min(m$P) * 100))

m3    <- fits_lca[[2]]
probs <- sapply(items, function(it) m3$probs[[it]][, 2])
round(m3$P, 3)
round(probs, 2)

df_lca$state <- match(m3$predclass, order(rowMeans(probs)))   # 1 = low ... 3 = high
table(df_lca$state, df_lca$wavename)


# ---- 4. Data for the latent Markov models ------------------------------------

df_lm <- make_lm_data(df_lca)

stopifnot(nrow(df_lm) == 46585)
table(table(df_lm$pidp))
colSums(df_lm[, emp_vars])   # 417 / 471 / 448
anyNA(df_lm[, c(emp_vars, "female", "age_c")])

# Transitions between modally assigned states (benchmark)
tt <- df_lm %>%
  group_by(pidp) %>%
  mutate(state_lag = lag(state), t_lag = lag(t)) %>%
  ungroup() %>%
  filter(!is.na(state_lag), t == t_lag + 1)
round(prop.table(table(from = tt$state_lag, to = tt$state), 1), 3)

tt2 <- df_lm %>%
  group_by(pidp) %>%
  mutate(state_lag = lag(state)) %>%
  ungroup() %>%
  filter(re_employed == 1, !is.na(state_lag))
table(tt2$state_lag)


# ---- 5. Helper functions -----------------------------------------------------

fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<.001", sprintf("%.3f", p)))
fmt_ms <- function(x) sprintf("%.1f (%.1f)", mean(x), sd(x))

state_order <- function(m) order(rowMeans(m$Psi[2, , ]))

mk_f <- function(covs) as.formula(paste("~ 1 |", paste(covs, collapse = " + ")))

fit_starts <- function(seeds, ...) {
  lapply(seeds, function(s) {
    set.seed(s)
    invisible(capture.output(m <- lmest(...)))
    m
  })
}

lr_test <- function(m0, m1) {
  stat <- 2 * (m1$lk - m0$lk)
  df   <- m1$np - m0$np
  c(LR = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE))
}

# Odds ratios; each logit compares moving to a state with staying in the origin
# state. Multiple-testing adjustment applies to the employment coefficients only.
or_table <- function(m, emp = emp_vars) {
  k    <- m$k
  labs <- if (k == 3) labs3 else paste0("S", 1:k)
  st   <- labs[rank(rowMeans(m$Psi[2, , ]))]
  Ga   <- m$Ga
  se   <- m$seGa
  out  <- do.call(rbind, lapply(1:k, function(r) {
    to <- st[setdiff(1:k, r)]
    do.call(rbind, lapply(1:(k - 1), function(cc) data.frame(
      From = st[r], To = to[cc], Covariate = dimnames(Ga)[[1]][-1],
      b = Ga[-1, cc, r], se = se[-1, cc, r])))
  }))
  out$OR <- exp(out$b)
  out$lo <- exp(out$b - 1.96 * out$se)
  out$hi <- exp(out$b + 1.96 * out$se)
  out$p  <- 2 * pnorm(-abs(out$b / out$se))
  e <- out$Covariate %in% emp
  out$p_holm <- NA
  out$p_fdr  <- NA
  out$p_holm[e] <- p.adjust(out$p[e], "holm")
  out$p_fdr[e]  <- p.adjust(out$p[e], "BH")
  out$From <- factor(out$From, labs)
  out$To   <- factor(out$To, labs)
  out[order(out$From, out$To, out$Covariate), ]
}

print_or <- function(tab) {
  print(cbind(tab[, c("From", "To", "Covariate")],
              round(tab[, c("OR", "lo", "hi", "p", "p_holm", "p_fdr")], 3)),
        row.names = FALSE)
}

pred_trans <- function(m, x) {
  k  <- m$k
  Ga <- m$Ga
  P  <- matrix(0, k, k)
  for (r in 1:k) {
    e <- c(1, exp(c(x %*% Ga[, , r])))
    p <- e / sum(e)
    P[r, r] <- p[1]
    P[r, setdiff(1:k, r)] <- p[-1]
  }
  o <- state_order(m)
  P <- P[o, o]
  dimnames(P) <- list(from = labs3, to = labs3)
  round(P, 3)
}

# Models A (no covariates), B (employment transitions), C (+ age, sex)
fit_models <- function(d, extra = NULL, seed_C = NULL) {
  fitsA <- fit_starts(1:6, responsesFormula = resp, index = idx, data = d,
                      k = 3, modBasic = 1, start = 1, tol = 1e-6, maxit = 2000)
  print(round(sapply(fitsA, `[[`, "lk"), 1))
  A <- fitsA[[which.max(sapply(fitsA, `[[`, "lk"))]]

  set.seed(123)
  B <- lmest(responsesFormula = resp, latentFormula = mk_f(c(emp_vars, extra)),
             index = idx, data = d, k = 3, start = 1, tol = 1e-6, maxit = 2000,
             out_se = TRUE)

  fC    <- mk_f(c(emp_vars, "female", "age_c", extra))
  fitsC <- fit_starts(1:5, responsesFormula = resp, latentFormula = fC,
                      index = idx, data = d, k = 3, start = 1, tol = 1e-6, maxit = 2000)
  lkC <- sapply(fitsC, `[[`, "lk")
  print(round(lkC, 1))
  set.seed(if (is.null(seed_C)) which.max(lkC) else seed_C)
  C <- lmest(responsesFormula = resp, latentFormula = fC,
             index = idx, data = d, k = 3, start = 1, tol = 1e-6, maxit = 2000,
             out_se = TRUE)
  list(A = A, B = B, C = C)
}

# Model C with stayed unemployed as the reference category
fit_ref <- function(d, extra = NULL) {
  set.seed(123)
  lmest(responsesFormula = resp,
        latentFormula = mk_f(c("stayed_emp", "job_loss", "re_employed", "female", "age_c", extra)),
        index = idx, data = d, k = 3, start = 1, tol = 1e-6, maxit = 2000, out_se = TRUE)
}

model_table <- function(ms, covariates = NULL) {
  tab <- data.frame(Model      = c("A", "B", "C"),
                    logLik     = sprintf("%.1f", sapply(ms, `[[`, "lk")),
                    Parameters = sapply(ms, `[[`, "np"),
                    AIC        = sprintf("%.1f", sapply(ms, `[[`, "aic")),
                    BIC        = sprintf("%.1f", sapply(ms, `[[`, "bic")))
  if (!is.null(covariates)) tab <- cbind(tab[1], Covariates = covariates, tab[-1])
  lr <- rbind(c(NA, NA, NA), lr_test(ms$A, ms$B), lr_test(ms$B, ms$C))
  tab$`LR (vs previous)` <- ifelse(is.na(lr[, 1]), "-", sprintf("%.1f", lr[, 1]))
  tab$df                 <- ifelse(is.na(lr[, 2]), "-", lr[, 2])
  tab$p                  <- ifelse(is.na(lr[, 3]), "-", fmt_p(lr[, 3]))
  tab
}


# ---- 6. Main analysis: balanced sample (all five waves) ----------------------

mods <- fit_models(df_lm, seed_C = 123)
mA <- mods$A
mB <- mods$B
mC <- mods$C

setNames(round(mA$piv[state_order(mA)], 3), labs3)   # initial distribution
PA <- mA$Pi[state_order(mA), state_order(mA), 2]
dimnames(PA) <- list(from = labs3, to = labs3)
round(PA, 3)

pred_trans(mB, c(1, 0, 0, 0))   # stayed employed
pred_trans(mB, c(1, 1, 0, 0))   # job loss
pred_trans(mB, c(1, 0, 1, 0))   # stayed unemployed
pred_trans(mB, c(1, 0, 0, 1))   # re-employment

# Timing check: covariates shifted one wave forward should fit worse
df_lead <- df_lm %>%
  arrange(pidp, t) %>%
  group_by(pidp) %>%
  mutate(across(all_of(emp_vars), ~ coalesce(lead(.), 0L), .names = "{.col}_lead")) %>%
  ungroup() %>%
  as.data.frame()

set.seed(123)
invisible(capture.output(
  mB_lead <- lmest(responsesFormula = resp,
                   latentFormula = ~ 1 | job_loss_lead + stayed_unemp_lead + re_employed_lead,
                   index = idx, data = df_lead, k = 3, start = 1, tol = 1e-6, maxit = 2000)))
c(model_B = mB$lk, lead_covariates = mB_lead$lk)

table2 <- model_table(mods, c("None", "Employment transitions",
                              "Employment transitions + age, sex"))
table2

tabB <- or_table(mB)
tabC <- or_table(mC)
print_or(tabC[tabC$Covariate %in% emp_vars, ])
print_or(tabC[!tabC$Covariate %in% emp_vars, ])

# Re-employment vs continued unemployment
emp_ref <- c("stayed_emp", "job_loss", "re_employed")
mC_ref  <- fit_ref(df_lm)
mC_ref$lk
tab_ref <- or_table(mC_ref, emp = emp_ref)
print_or(tab_ref[tab_ref$Covariate == "re_employed", ])

re_ref <- tab_ref[tab_ref$Covariate == "re_employed", ]
table3_ref <- data.frame(
  Transition    = paste(re_ref$From, "to", re_ref$To),
  `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", re_ref$OR, re_ref$lo, re_ref$hi),
  p             = fmt_p(re_ref$p),
  `p (Holm)`    = fmt_p(re_ref$p_holm),
  `p (BH)`      = fmt_p(re_ref$p_fdr),
  check.names = FALSE)
table3_ref

# Employment effects before and after adjustment for age and sex
cmp_adj <- merge(tabB[tabB$Covariate %in% emp_vars, c("From", "To", "Covariate", "OR", "p_holm")],
                 tabC[tabC$Covariate %in% emp_vars, c("From", "To", "Covariate", "OR", "lo", "hi", "p_holm", "p_fdr")],
                 by = c("From", "To", "Covariate"), suffixes = c("_B", "_C"))
cmp_adj <- cmp_adj[order(cmp_adj$From, cmp_adj$To, cmp_adj$Covariate), ]
print(cbind(cmp_adj[, 1:3], round(cmp_adj[, -(1:3)], 3)), row.names = FALSE)


# ---- 7. Sensitivity analysis I: four latent states ---------------------------

fC <- mk_f(c(emp_vars, "female", "age_c"))
fits4 <- fit_starts(1:5, responsesFormula = resp, latentFormula = fC,
                    index = idx, data = df_lm, k = 4, start = 1, tol = 1e-6, maxit = 2000)
lk4 <- sapply(fits4, `[[`, "lk")
round(lk4, 1)

set.seed(which.max(lk4))
mC4 <- lmest(responsesFormula = resp, latentFormula = fC,
             index = idx, data = df_lm, k = 4, start = 1, tol = 1e-6, maxit = 2000,
             out_se = TRUE)
c(best_start = max(lk4), with_se = mC4$lk)
c(K3_AIC = mC$aic, K4_AIC = mC4$aic, K3_BIC = mC$bic, K4_BIC = mC4$bic)

o4    <- state_order(mC4)
labs4 <- paste0("S", 1:4)
prof4 <- mC4$Psi[2, o4, ]
dimnames(prof4) <- list(state = labs4, item = items)
round(prof4, 2)
setNames(round(colMeans(mC4$Piv)[o4], 3), labs4)

tab4 <- or_table(mC4)
print_or(tab4[tab4$Covariate %in% emp_vars, ])


# ---- 8. Sensitivity analysis II: unbalanced panel (>= 3 waves) ---------------

df_u <- make_lm_data(prepare_sample(min_waves = 3))

cat("Individuals:", n_distinct(df_u$pidp), "\n")
cat("Person-waves:", nrow(df_u), "\n")
table(waves_per_person = table(df_u$pidp))
colSums(df_u[, c(emp_vars, "after_gap")])
round(100 * mean(df_u$emp_status == "unemployed"), 1)
stopifnot(!anyNA(df_u[, c(emp_vars, "after_gap", "female", "age_c")]))

mods_u <- fit_models(df_u, extra = "after_gap")
uA <- mods_u$A
uB <- mods_u$B
uC <- mods_u$C

tableE2 <- model_table(mods_u)
tableE2

tab_u <- or_table(uC)
key   <- function(t) t[t$Covariate %in% emp_vars, c("From", "To", "Covariate", "OR", "lo", "hi")]
cmp_u <- merge(key(tabC), key(tab_u), by = c("From", "To", "Covariate"),
               suffixes = c("_balanced", "_unbalanced"))
cmp_u <- cmp_u[order(cmp_u$From, cmp_u$To, cmp_u$Covariate), ]
print(cbind(cmp_u[, 1:3], round(cmp_u[, -(1:3)], 2)), row.names = FALSE)

uC_ref   <- fit_ref(df_u, extra = "after_gap")
tab_ref_u <- or_table(uC_ref, emp = emp_ref)
print_or(tab_ref_u[tab_ref_u$Covariate == "re_employed", ])


# ---- 9. Sensitivity analyses III-V -------------------------------------------
# III: lead (placebo) test, IV: duration of unemployment, V: employment-sensitive
# items removed

fit_C <- function(d, covs, resp_f = resp, seeds = 1:5) {
  f    <- mk_f(covs)
  fits <- fit_starts(seeds, responsesFormula = resp_f, latentFormula = f,
                     index = idx, data = d, k = 3, start = 1, tol = 1e-6, maxit = 2000)
  lk <- sapply(fits, `[[`, "lk")
  print(round(lk, 1))
  set.seed(which.max(lk))
  lmest(responsesFormula = resp_f, latentFormula = f, index = idx, data = d,
        k = 3, start = 1, tol = 1e-6, maxit = 2000, out_se = TRUE)
}

fmt_or <- function(t) sprintf("%.2f (%.2f-%.2f)", t$OR, t$lo, t$hi)

# III. Lead test (balanced sample): do next wave's employment transitions predict
# distress transitions? Joint test of the 18 lead coefficients against Model C
# with the same lead_missing indicator but no lead dummies.
lead_vars <- paste0(emp_vars, "_lead")
colSums(df_lm[, c(lead_vars, "lead_missing")])

mC_lm   <- fit_C(df_lm, c(emp_vars, "female", "age_c", "lead_missing"))
mC_lead <- fit_C(df_lm, c(emp_vars, "female", "age_c", "lead_missing", lead_vars))

lr_test(mC_lm, mC_lead)

tab_lead <- or_table(mC_lead, emp = lead_vars)
print_or(tab_lead[tab_lead$Covariate %in% lead_vars, ])

lead_rows <- tab_lead[tab_lead$Covariate %in% lead_vars, ]
tableF1 <- data.frame(
  Transition    = paste(lead_rows$From, "to", lead_rows$To),
  Covariate     = sub("_lead$", "", as.character(lead_rows$Covariate)),
  `OR (95% CI)` = fmt_or(lead_rows),
  p             = fmt_p(lead_rows$p),
  `p (Holm)`    = fmt_p(lead_rows$p_holm),
  check.names = FALSE)
tableF1

# Current-wave effects with and without the lead terms
cur <- function(t) t[t$Covariate %in% emp_vars, c("From", "To", "Covariate", "OR", "lo", "hi")]
cmp_lead <- merge(cur(tabC), cur(tab_lead), by = c("From", "To", "Covariate"),
                  suffixes = c("_main", "_with_lead"))
cmp_lead <- cmp_lead[order(cmp_lead$From, cmp_lead$To, cmp_lead$Covariate), ]
print(cbind(cmp_lead[, 1:3], round(cmp_lead[, -(1:3)], 2)), row.names = FALSE)

# IV. Duration of unemployment (unbalanced sample). Spells that began before
# Wave 6 cannot be dated, so "2w" may include longer spells (left-censoring).
emp_dur <- c("job_loss", "stayed_unemp_2w", "stayed_unemp_3w", "re_employed")
colSums(df_u[, c(emp_dur, "after_gap")])

uC_dur  <- fit_C(df_u, c(emp_dur, "female", "age_c", "after_gap"))
tab_dur <- or_table(uC_dur, emp = emp_dur)
print_or(tab_dur[tab_dur$Covariate %in% c("stayed_unemp_2w", "stayed_unemp_3w"), ])

# Third or later wave vs. second wave of unemployment (same model, reference changed)
emp_dur_ref <- c("stayed_emp", "job_loss", "stayed_unemp_3w", "re_employed")
uC_dur_ref  <- fit_C(df_u, c(emp_dur_ref, "female", "age_c", "after_gap"))
c(main = uC_dur$lk, reference_changed = uC_dur_ref$lk)   # should be equal
tab_dur_ref <- or_table(uC_dur_ref, emp = emp_dur_ref)
print_or(tab_dur_ref[tab_dur_ref$Covariate == "stayed_unemp_3w", ])

dur_rows <- tab_dur[tab_dur$Covariate %in% c("stayed_unemp_2w", "stayed_unemp_3w"), ]
tableF2 <- data.frame(
  Transition    = paste(dur_rows$From, "to", dur_rows$To),
  Covariate     = ifelse(dur_rows$Covariate == "stayed_unemp_2w",
                         "2nd wave unemployed", "3rd or later wave unemployed"),
  `OR (95% CI)` = fmt_or(dur_rows),
  p             = fmt_p(dur_rows$p),
  `p (Holm)`    = fmt_p(dur_rows$p_holm),
  check.names = FALSE)
tableF2 <- tableF2[order(dur_rows$From, dur_rows$To, dur_rows$Covariate), ]
tableF2

# V. Items most tied to employment removed (c useful role, g enjoyment)
items_drop <- c("scghqc", "scghqg")
resp_drop  <- as.formula(paste(paste(setdiff(items, items_drop), collapse = " + "), "~ NULL"))

mC_drop <- fit_C(df_lm, c(emp_vars, "female", "age_c"), resp_f = resp_drop)
round(mC_drop$Psi[2, state_order(mC_drop), ], 2)   # profiles: same three states?
round(setNames(colMeans(mC_drop$Piv)[state_order(mC_drop)], labs3), 3)

tab_drop <- or_table(mC_drop)
cmp_drop <- merge(cur(tabC), cur(tab_drop), by = c("From", "To", "Covariate"),
                  suffixes = c("_main", "_10items"))
cmp_drop <- cmp_drop[order(cmp_drop$From, cmp_drop$To, cmp_drop$Covariate), ]
print(cbind(cmp_drop[, 1:3], round(cmp_drop[, -(1:3)], 2)), row.names = FALSE)

tableF3 <- cmp_drop
tableF3[, -(1:3)] <- round(tableF3[, -(1:3)], 2)


# ---- 10. Tables for the paper ------------------------------------------------

describe <- function(d) {
  d %>% summarise(
    `Individuals`                   = format(n_distinct(pidp), big.mark = ","),
    `Person-waves`                  = format(n(), big.mark = ","),
    `Unemployed (%)`                = sprintf("%.1f", 100 * mean(emp_status == "unemployed")),
    `Age, mean (SD)`                = fmt_ms(age_dv),
    `Women (%)`                     = sprintf("%.1f", 100 * mean(female)),
    `GHQ caseness score, mean (SD)` = fmt_ms(scghq2_dv),
    `GHQ caseness score >= 4 (%)`   = sprintf("%.1f", 100 * mean(scghq2_dv >= 4)),
    `Job loss (n)`                  = as.character(sum(job_loss)),
    `Stayed unemployed (n)`         = as.character(sum(stayed_unemp)),
    `Re-employment (n)`             = as.character(sum(re_employed)))
}

t1 <- bind_rows(
  df_lm %>% group_by(wavename) %>% describe() %>% ungroup() %>%
    mutate(Wave = paste("Wave", wavename)) %>% dplyr::select(-wavename),
  describe(df_lm) %>% mutate(Wave = "Total"))
t1 <- t1[, c("Wave", setdiff(names(t1), "Wave"))]
t1[t1$Wave == "Wave 6", c("Job loss (n)", "Stayed unemployed (n)", "Re-employment (n)")] <- "-"
table1 <- data.frame(Characteristic = names(t1)[-1], t(t1[, -1]), check.names = FALSE)
names(table1)[-1] <- t1$Wave
table1

cov_labels <- c(job_loss = "Job loss", stayed_unemp = "Stayed unemployed",
                re_employed = "Re-employment", female = "Women", age_c = "Age (per decade)")
table3 <- data.frame(
  Transition    = paste(tabC$From, "to", tabC$To),
  Covariate     = cov_labels[as.character(tabC$Covariate)],
  `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", tabC$OR, tabC$lo, tabC$hi),
  p             = fmt_p(tabC$p),
  `p (Holm)`    = fmt_p(tabC$p_holm),
  `p (BH)`      = fmt_p(tabC$p_fdr),
  check.names = FALSE)
table3 <- table3[order(factor(tabC$From, labs3), factor(tabC$To, labs3),
                       factor(as.character(tabC$Covariate), names(cov_labels))), ]
table3

# Appendix tables: balanced vs unbalanced sample
tableE1 <- data.frame(
  Characteristic = c("Individuals", "Person-waves", "Unemployed (%)",
                     "Job loss (n)", "Stayed unemployed (n)", "Re-employment (n)"),
  Balanced   = c(format(n_distinct(df_lm$pidp), big.mark = ","),
                 format(nrow(df_lm), big.mark = ","),
                 sprintf("%.1f", 100 * mean(df_lm$emp_status == "unemployed")),
                 colSums(df_lm[, emp_vars])),
  Unbalanced = c(format(n_distinct(df_u$pidp), big.mark = ","),
                 format(nrow(df_u), big.mark = ","),
                 sprintf("%.1f", 100 * mean(df_u$emp_status == "unemployed")),
                 colSums(df_u[, emp_vars])))
tableE1

tableE3 <- cmp_u
tableE3[, -(1:3)] <- round(tableE3[, -(1:3)], 2)
tableE3


# ---- 11. Figures -------------------------------------------------------------

state_cols <- c(Low = "deepskyblue", Moderate = "blue", High = "black")
theme_paper <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        axis.text  = element_text(colour = "grey30"),
        legend.position = "top")

# Figure 1: probability of a distress response by item and latent state
item_labels <- c("a Concentration", "b Sleep loss", "c Useful role", "d Decisions",
                 "e Under strain", "f Overcoming difficulties", "g Enjoyment",
                 "h Facing problems", "i Unhappy", "j Confidence", "k Worthless",
                 "l Happiness")
prof <- mC$Psi[2, state_order(mC), ]
round(prof, 2)
fig1_dat <- data.frame(
  state = factor(rep(labs3, times = 12), labs3),
  item  = factor(rep(item_labels, each = 3), item_labels),
  prob  = as.vector(prof))

fig1 <- ggplot(fig1_dat, aes(item, prob, colour = state, group = state, shape = state)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  geom_text(data = subset(fig1_dat, item == "l Happiness"),
            aes(label = state), hjust = -0.25, size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = state_cols, name = "Latent state") +
  scale_shape_manual(values = c(Low = 16, Moderate = 17, High = 15), name = "Latent state") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.4))) +
  labs(x = "GHQ-12 item", y = "Probability of a distress response") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
fig1

# Figure 2: average predicted transition probabilities by employment transition
avg_trans <- function(m, data, scenario) {
  Ga  <- m$Ga
  k   <- m$k
  cov <- dimnames(Ga)[[1]][-1]
  X   <- cbind(1, as.matrix(data[, cov]))
  colnames(X) <- c("(Intercept)", cov)
  X[, emp_vars] <- 0
  if (scenario %in% emp_vars) X[, scenario] <- 1
  P <- matrix(0, k, k)
  for (r in 1:k) {
    e <- cbind(1, exp(X %*% Ga[, , r]))
    p <- colMeans(e / rowSums(e))
    P[r, r] <- p[1]
    P[r, setdiff(1:k, r)] <- p[-1]
  }
  o <- state_order(m)
  P <- P[o, o]
  dimnames(P) <- list(from = labs3, to = labs3)
  P
}

scenarios <- c(stayed_employed = "Stayed employed", job_loss = "Job loss",
               stayed_unemp = "Stayed unemployed", re_employed = "Re-employment")
panels <- c("Low to Moderate or High\n(deterioration, H1)",
            "High to Low\n(recovery, H2)",
            "Moderate to Low\n(recovery, H3)")

fig2_dat <- do.call(rbind, lapply(names(scenarios), function(s) {
  P <- avg_trans(mC, df_lm, s)
  data.frame(group    = scenarios[[s]],
             quantity = panels,
             prob     = c(P["Low", "Moderate"] + P["Low", "High"],
                          P["High", "Low"],
                          P["Moderate", "Low"]))
}))
fig2_dat$group    <- factor(fig2_dat$group, rev(scenarios))
fig2_dat$quantity <- factor(fig2_dat$quantity, panels)
fig2_dat

fig2 <- ggplot(fig2_dat, aes(prob, group)) +
  geom_segment(aes(x = 0, xend = prob, yend = group), colour = "grey", linewidth = 0.6) +
  geom_point(colour = "blue", size = 2.8) +
  geom_text(aes(label = sprintf("%.2f", prob)), hjust = -0.4, size = 3.4, colour = "black") +
  facet_wrap(~ quantity, nrow = 1) +
  scale_x_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "Predicted transition probability", y = NULL) +
  theme_paper +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
fig2


# ---- 12. Save ----------------------------------------------------------------

out <- function(f) file.path(project, f)

write.csv(table1,   out("Table1_sample.csv"),                 row.names = FALSE)
write.csv(table2,   out("Table2_model_comparison.csv"),       row.names = FALSE)
write.csv(table3,   out("Table3_odds_ratios.csv"),            row.names = FALSE)
write.csv(table3_ref, out("Table3_reemployment_vs_stayed_unemployed.csv"), row.names = FALSE)
write.csv(tabC,     out("table_modelC_K3.csv"),               row.names = FALSE)
write.csv(tab4,     out("table_modelC_K4.csv"),               row.names = FALSE)
write.csv(tableE1,  out("TableE1_sample_comparison.csv"),     row.names = FALSE)
write.csv(tableE2, out("TableE2_model_comparison_unbalanced.csv"), row.names = FALSE)
write.csv(tableE3,  out("TableE3_balanced_vs_unbalanced.csv"), row.names = FALSE)
write.csv(tab_u,    out("table_modelC_unbalanced.csv"),       row.names = FALSE)
write.csv(tableF1,  out("TableF1_lead_test.csv"),             row.names = FALSE)
write.csv(cmp_lead, out("TableF1b_current_effects_with_lead.csv"), row.names = FALSE)
write.csv(tableF2,  out("TableF2_unemployment_duration.csv"), row.names = FALSE)
write.csv(tableF3,  out("TableF3_ten_items.csv"),             row.names = FALSE)

ggsave(out("Figure1_item_profiles.png"), fig1, width = 7.5, height = 4.5, dpi = 300)
ggsave(out("Figure1_item_profiles.pdf"), fig1, width = 7.5, height = 4.5)
ggsave(out("Figure2_transitions.png"),   fig2, width = 8,   height = 3.2, dpi = 300)
ggsave(out("Figure2_transitions.pdf"),   fig2, width = 8,   height = 3.2)

saveRDS(list(mA = mA, mB = mB, mC = mC, mC4 = mC4), out("lm_models.rds"))
saveRDS(list(uA = uA, uB = uB, uC = uC, uC_ref = uC_ref), out("lm_models_unbalanced.rds"))
saveRDS(list(mC_lm = mC_lm, mC_lead = mC_lead, uC_dur = uC_dur,
             uC_dur_ref = uC_dur_ref, mC_drop = mC_drop), out("lm_models_extra.rds"))

