# =============================================================================
# swan_rule_diagnostics.R   (v3: full covariate set, ordinal HOTFLAS)
#
# Outcome-blind support diagnostics for dynamic MHT-initiation rules in SWAN:
#   (1) adherence table, (2) positivity checks, (3) cumulative-weight
#   diagnostics, (4) pooled event counts, (5) model-size (events-per-parameter)
#   checks for the g models, (6) data-adaptive choice of K and delta.
#
# Nothing here fits an outcome (Q) regression or tabulates events by rule.
#
# ---- Data layout expected (wide, one row per woman) -------------------------
#   SWANID, AGE, RACE.* dummies, DEGREE,                       (baseline)
#   <stem>0 for each time-varying stem, MHT0, HTN0,            (visit 0)
#   <stem>v, D_TVCOV<v>, D_MHT<v>, MHT<v>, C_HTN<v>, HTN<v>    (visits v = 1..K)
# Time-varying stems in your file: DIABETE, PHYSACT, HEIGHT, WEIGHT, STATUS,
# SMOKERE, HOTFLAS, DRNKBEE, DEPRESSION, ANXIETY.
#
# Codings (from the SWAN codebook you supplied):
#   STATUS : 0 post by BSO | 1 natural post | 2 late peri | 3 early peri | 4 pre |
#            5 pregnant/breastfeeding | 6 unknown due to HT use | 7 unknown due to hysterectomy
#   HOTFLAS: 0 not at all | 1 1-5 days | 2 6-8 days | 3 9-13 days | 4 every day
#            (past two weeks). ORDINAL, not binary: set hot_codes explicitly.
#   MHT, HTN, D_MHT, D_TVCOV, C_HTN are 0/1.
#   C_HTN = 1 -> censored at that visit. In your file, after censoring the covariates
#   are carried forward, D_* flags are 1 and HTN is 0 (not NA); the script handles this
#   because events are only counted while uncensored.
#   D_MHT = 1 -> MHT variable missing at that visit.
#
# ---- Node ordering ----------------------------------------------------------
#   W = baseline covariates, then for v = 1..K:
#   A_v = MHT_v (use since last visit), C_v, Y_v = HTN_v, L_v = time-varying covariates
#   so A_v and C_v are modeled on covariates from visit v-1 (lagged one visit).
#
# ---- Rules (all "initiate and continue") ------------------------------------
#   never    : A_v = 0 for all v
#   Stage-triggered (needs status_trigger): tau_S = first visit v >= 0 with
#   STATUS_v in status_trigger
#     early    : A_v = 1 iff tau_S <= v - 1
#     late_dK  : A_v = 1 iff tau_S <= v - 1 - K     (K visits after the trigger)
#   Symptom-triggered (hot_rule = TRUE): tau_H = first visit v >= 0 with HOTFLAS_v in hot_codes
#     hot_early: A_v = 1 iff tau_H <= v - 1
#     hot_dK   : A_v = 1 iff tau_H <= v - 1 - K     (optional, via deltas_hot)
#   Because covariates at v-1 are observed before A_v, "as soon as symptoms are
#   reported" means initiating in the interval after the visit at which they were
#   first reported. Once started, MHT is continued; women who stop are artificially censored.
#
#   Suggested codes: status_trigger = c(2, 1) for "late perimenopause or later";
#   hot_codes = 1:4 for "any hot flashes in the past two weeks" (hot_codes = 3:4 for
#   "frequent": 9+ days). Run hot_codes both ways as a sensitivity analysis.
#
# ---- Choices you may want to change -----------------------------------------
#   * Missing MHT (D_MHT = 1, or MHT NA) and missing HTN are treated as censoring.
#   * g_C is modeled on history through visit v-1 (not on A_v, which can be missing).
#     ltmle conditions the censoring node on A_v; this is a simplification.
#   * g models are fit on all women at risk (and uncensored, for g_A) at the node,
#     with A_{v-1} as a covariate.
#   * IMPORTANT: positivity and weight diagnostics are only meaningful for the
#     adjustment set you will use in the real analysis. Run with the full set.
#   * Default learner is main-effects logistic regression (glm.fit). Use fit_fun
#     (treatment) and fit_fun_C (censoring) to swap in a penalized learner.
#
# ---- Changes from v2 --------------------------------------------------------
#   * Full covariate set in the g models (see `covs` / swan_covariates()):
#     numeric stems are lagged, last-observation-carried-forward, standardized;
#     categorical stems (STATUS, DEGREE) become factors with an explicit "NA" level;
#     lagged D_TVCOV and a "any covariate missing" flag are included.
#   * `rhs` is gone; use `covs` instead.
#   * hot_codes: HOTFLAS is ordinal, so you must say which values count as symptoms.
#     (Auto-set to 1 only if HOTFLAS is binary; otherwise the function stops.)
#   * model_info table + warning when a g model has fewer than epv_min events (or
#     non-events) per fitted parameter at some visit.
#   * fit_fun_C: separate learner for the censoring model (few women are censored
#     at each visit, so it overfits first). fit_glmnet() provided (needs glmnet).
#   * exclude_status0: drop women with certain baseline STATUS codes (e.g. c(0, 7, 5)).
# =============================================================================

# ---------- small helpers ----------------------------------------------------

get_mat <- function(dat, prefix, K) {
  M <- matrix(NA_real_, nrow(dat), K + 1)
  for (v in 0:K) {
    nm <- paste0(prefix, v)
    if (nm %in% names(dat)) M[, v + 1] <- as.numeric(dat[[nm]])
  }
  M
}

# first visit in v_seq at which logical matrix M (columns = visits 0..K) is TRUE
first_time <- function(M, v_seq) {
  out <- rep(Inf, nrow(M))
  for (v in rev(v_seq)) out[which(M[, v + 1])] <- v
  out
}

ess_fun <- function(w) if (length(w) == 0) NA_real_ else sum(w)^2 / sum(w^2)
qfun <- function(x, p) if (length(x) == 0) NA_real_ else as.numeric(quantile(x, p, na.rm = TRUE))
mfun <- function(x, f) if (length(x) == 0) NA_real_ else f(x)

auc_fun <- function(y, p) {
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Default learner: logistic regression on a numeric design matrix.
# Aliased columns (e.g., empty factor levels, collinear terms) get coefficient 0.
fit_glm <- function(y, X_fit, X_new) {
  if (length(unique(y)) < 2) return(rep(mean(y), nrow(X_new)))
  f <- suppressWarnings(glm.fit(X_fit, y, family = binomial()))
  b <- f$coefficients; b[is.na(b)] <- 0
  as.numeric(plogis(X_new %*% b))
}

# Penalized alternative (ridge by default; alpha = 1 gives lasso). Needs the glmnet
# package. Falls back to the marginal mean when there are too few events/non-events.
fit_glmnet <- function(y, X_fit, X_new, alpha = 0) {
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Install glmnet first: install.packages('glmnet')")
  nmin <- min(sum(y == 1), sum(y == 0))
  if (nmin < 3) return(rep(mean(y), nrow(X_new)))
  tryCatch({
    cv <- suppressWarnings(glmnet::cv.glmnet(X_fit[, -1, drop = FALSE], y, family = "binomial",
                                             alpha = alpha, nfolds = min(5, nmin)))
    as.numeric(predict(cv, X_new[, -1, drop = FALSE], s = "lambda.min", type = "response"))
  }, error = function(e) rep(mean(y), nrow(X_new)))
}

cv_auc <- function(y, X, fit_fun, folds = 5) {
  n <- length(y)
  if (n < 2 * folds || length(unique(y)) < 2) return(NA_real_)
  id <- sample(rep(seq_len(folds), length.out = n)); p <- numeric(n)
  for (k in seq_len(folds)) {
    tr <- id != k
    p[!tr] <- fit_fun(y[tr], X[tr, , drop = FALSE], X[!tr, , drop = FALSE])
  }
  auc_fun(y, p)
}

# Default covariate set, matched to the column names in the SWAN extract.
# Override by passing covs = list(W_num=, W_fac=, L_num=, L_fac=) to run_diagnostics().
#   W_num: baseline numeric columns        W_fac: baseline categorical columns
#   L_num: time-varying numeric stems      L_fac: time-varying categorical stems
swan_covariates <- function(dat) {
  nm <- names(dat)
  list(
    W_num = intersect(c("AGE", grep("^RACE", nm, value = TRUE)), nm),
    W_fac = intersect("DEGREE", nm),
    L_num = Filter(function(s) paste0(s, "0") %in% nm,
                   c("DIABETE", "PHYSACT", "HEIGHT", "WEIGHT", "SMOKERE",
                     "DRNKBEE", "DEPRESSION", "ANXIETY", "HOTFLAS")),
    L_fac = intersect("STATUS", sub("[0-9]+$", "", nm)))
}

# ---------- main function ----------------------------------------------------

run_diagnostics <- function(dat, K = 10,
                            status_trigger = NULL, deltas = c(1, 2),
                            hot_rule = TRUE, hot_codes = NULL, deltas_hot = NULL,
                            exclude_status0 = NULL,
                            covs = NULL,
                            bounds = c(0.01, 0.05),
                            fit_fun = fit_glm, fit_fun_C = fit_fun,
                            cv_auc_folds = 5, epv_min = 10,
                            outdir = NULL, seed = 1) {
  if (is.null(status_trigger) && !hot_rule)
    stop("Nothing to evaluate: set status_trigger and/or hot_rule = TRUE.")
  set.seed(seed)
  
  # new-user, hypertension-free at baseline (optionally also drop some baseline STATUS codes)
  dat <- dat[which(dat$MHT0 %in% 0 & dat$HTN0 %in% 0), , drop = FALSE]
  if (!is.null(exclude_status0)) dat <- dat[!(dat$STATUS0 %in% exclude_status0), , drop = FALSE]
  n <- nrow(dat)
  message("Analytic sample (MHT-naive, HTN-free at baseline",
          if (!is.null(exclude_status0)) paste0(", STATUS0 not in {", paste(exclude_status0, collapse = ","), "}"),
          "): n = ", n)
  
  # ---- covariate specification ----
  if (is.null(covs)) covs <- swan_covariates(dat)
  covs <- modifyList(list(W_num = character(0), W_fac = character(0),
                          L_num = character(0), L_fac = character(0)), covs)
  bad <- setdiff(c(covs$W_num, covs$W_fac), names(dat))
  if (length(bad)) stop("Baseline covariates not found in data: ", paste(bad, collapse = ", "))
  for (s in c(covs$L_num, covs$L_fac))
    if (!any(paste0(s, 0:K) %in% names(dat))) stop("No columns found for time-varying stem: ", s)
  
  S  <- get_mat(dat, "STATUS", K);  Hf <- get_mat(dat, "HOTFLAS", K)
  A  <- get_mat(dat, "MHT", K);     Dm <- get_mat(dat, "D_MHT", K)
  Cc <- get_mat(dat, "C_HTN", K);   Y  <- get_mat(dat, "HTN", K)
  A[, 1] <- 0
  
  # ---- trigger definitions ----
  if (!is.null(status_trigger)) {
    unseen <- setdiff(status_trigger, unique(S[!is.na(S)]))
    if (length(unseen)) warning("status_trigger codes never observed: ", paste(unseen, collapse = ", "))
  }
  if (hot_rule) {
    hv <- sort(unique(Hf[!is.na(Hf)]))
    if (is.null(hot_codes)) {
      if (all(hv %in% c(0, 1))) hot_codes <- 1
      else stop("HOTFLAS takes values {", paste(hv, collapse = ","), "}, so 'HOTFLAS == 1' is ambiguous. ",
                "Set hot_codes, e.g. 1:4 for any hot flashes or 3:4 for 9+ days.")
    }
    if (!any(hv %in% hot_codes)) stop("hot_codes never observed in the data.")
    message("Symptom trigger: HOTFLAS in {", paste(hot_codes, collapse = ","), "}")
  }
  
  # censoring indicator at visit v: LTFU, missing MHT, or missing HTN
  cstar <- matrix(0, n, K + 1)
  for (v in 1:K) {
    cs <- (Cc[, v + 1] %in% 1) | (Dm[, v + 1] %in% 1) |
      is.na(A[, v + 1]) | is.na(Y[, v + 1])
    cstar[, v + 1] <- as.numeric(cs)
  }
  cens_time <- first_time(cstar == 1, 1:K)
  
  ev <- matrix(FALSE, n, K + 1)
  for (v in 1:K) ev[, v + 1] <- (cstar[, v + 1] == 0) & (Y[, v + 1] %in% 1)
  ev_time <- first_time(ev, 1:K)
  ev_time[cens_time < ev_time] <- Inf   # event after censoring is unobserved
  
  # node-v indicators: active = event-free and uncensored through v-1;
  # unc = additionally uncensored at v
  act <- sapply(1:K, function(v) ev_time >= v & cens_time >= v)
  unc <- sapply(1:K, function(v) cens_time > v)
  
  # trigger visits (first visit at which the trigger is observed, while still
  # event-free and uncensored)
  lim <- pmin(ev_time, cens_time - 1, K)
  first_trigger <- function(trigM) {
    tau <- rep(Inf, n)
    for (v in rev(0:(K - 1))) tau[which(trigM[, v + 1] & (v <= lim))] <- v
    tau
  }
  tau_S <- tau_H <- rep(Inf, n)
  if (!is.null(status_trigger)) {
    trigS <- matrix(FALSE, n, K + 1)
    for (v in 0:K) trigS[, v + 1] <- S[, v + 1] %in% status_trigger
    tau_S <- first_trigger(trigS)
  }
  if (hot_rule) {
    trigH <- matrix(FALSE, n, K + 1)
    for (v in 0:K) trigH[, v + 1] <- Hf[, v + 1] %in% hot_codes
    tau_H <- first_trigger(trigH)
  }
  
  # ---- rules: each is list(tau = trigger visit per woman, delta = delay) ----
  rules <- list(never = list(tau = rep(Inf, n), delta = 0))
  if (!is.null(status_trigger)) {
    rules$early <- list(tau = tau_S, delta = 0)
    for (dl in deltas) rules[[paste0("late_d", dl)]] <- list(tau = tau_S, delta = dl)
  }
  if (hot_rule) {
    rules$hot_early <- list(tau = tau_H, delta = 0)
    for (dl in deltas_hot) rules[[paste0("hot_d", dl)]] <- list(tau = tau_H, delta = dl)
  }
  
  make_d <- function(tau, delta) {
    M <- matrix(0, n, K)
    for (v in 1:K) M[, v] <- as.numeric(tau <= v - 1 - delta)
    M
  }
  
  rule_obj <- lapply(rules, function(r) {
    d <- make_d(r$tau, r$delta)
    dev_time <- rep(Inf, n)
    for (v in 1:K) {
      idx <- which(act[, v] & unc[, v] & (A[, v + 1] != d[, v]) & is.infinite(dev_time))
      dev_time[idx] <- v
    }
    stop_time <- pmin(cens_time, dev_time)
    # I[i, v] = followed the rule and stayed uncensored through v
    # (weights carry forward after the absorbing event)
    I <- sapply(1:K, function(v) stop_time > pmin(v, ev_time))
    list(d = d, I = I)
  })
  
  # ---- design matrices for g_A and g_C (rule-independent) ----
  fac <- function(x, lev)
    factor(ifelse(is.na(x), "NA", as.character(x)), levels = c(as.character(lev), "NA"))
  
  Wdf <- data.frame(row.names = seq_len(n))
  for (s in covs$W_num) {
    x <- as.numeric(dat[[s]])
    if (anyNA(x)) { warning("Baseline covariate ", s, " has NAs; mean-imputed."); x[is.na(x)] <- mean(x, na.rm = TRUE) }
    sdv <- sd(x); if (!is.finite(sdv) || sdv == 0) sdv <- 1
    Wdf[[paste0("W_", s)]] <- (x - mean(x)) / sdv
  }
  for (s in covs$W_fac) {
    x <- dat[[s]]
    Wdf[[paste0("W_", s)]] <- fac(x, sort(unique(x[!is.na(x)])))
  }
  
  Znum <- list(); TVNA <- matrix(FALSE, n, K + 1)
  for (s in covs$L_num) {
    M <- get_mat(dat, s, K); TVNA <- TVNA | is.na(M)
    Fm <- M                                              # last observation carried forward
    for (u in 2:(K + 1)) { m <- is.na(Fm[, u]); Fm[m, u] <- Fm[m, u - 1] }
    mu <- mean(Fm, na.rm = TRUE); sdv <- sd(as.vector(Fm), na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    Z <- (Fm - mu) / sdv; Z[is.na(Z)] <- 0               # still-missing (no prior value) -> mean
    Znum[[s]] <- Z
  }
  Mfac <- list(); Lev <- list()
  for (s in covs$L_fac) {
    Mfac[[s]] <- get_mat(dat, s, K)
    Lev[[s]] <- sort(unique(Mfac[[s]][!is.na(Mfac[[s]])]))
  }
  use_tvna <- any(TVNA)
  use_dtv <- "D_TVCOV1" %in% names(dat)
  Dt <- get_mat(dat, "D_TVCOV", K); Dt[is.na(Dt)] <- 0
  
  design <- function(v) {                                # column v of each matrix = visit v-1
    Xdf <- Wdf
    for (s in names(Znum)) Xdf[[paste0(s, "_prev")]] <- Znum[[s]][, v]
    for (s in names(Mfac)) Xdf[[paste0(s, "_prev")]] <- fac(Mfac[[s]][, v], Lev[[s]])
    Xdf$MHT_prev <- ifelse(is.na(A[, v]), 0, A[, v])
    if (!is.null(status_trigger)) Xdf$ever_trig_prev <- as.numeric(tau_S <= v - 1)
    if (hot_rule)                 Xdf$ever_hot_prev  <- as.numeric(tau_H <= v - 1)
    if (use_tvna) Xdf$tv_na_prev <- as.numeric(TVNA[, v])
    if (use_dtv)  Xdf$d_tvcov_prev <- Dt[, v]
    model.matrix(~ ., data = Xdf)
  }
  
  pA <- matrix(1, n, K); pC <- matrix(1, n, K); auc <- rep(NA_real_, K)
  minfo <- vector("list", K)
  for (v in 1:K) {
    a <- which(act[, v]); u <- a[unc[a, v]]
    if (length(a) < 20 || length(u) < 20) {
      warning("Visit ", v, ": fewer than 20 women at risk; g set to 1. Interpret nothing past here.")
      minfo[[v]] <- data.frame(visit = v, n_active = length(a), n_unc = length(u), n_params = NA,
                               n_A1 = NA, n_A0 = NA, epv_A = NA, n_cens = NA, epv_C = NA)
      next
    }
    X <- design(v)
    p <- sum(apply(X[a, -1, drop = FALSE], 2, function(z) length(unique(z)) > 1)) + 1
    yC <- 1 - cstar[a, v + 1]
    pC[a, v] <- fit_fun_C(yC, X[a, , drop = FALSE], X[a, , drop = FALSE])
    yA <- A[u, v + 1]
    pA[a, v] <- fit_fun(yA, X[u, , drop = FALSE], X[a, , drop = FALSE])
    if (cv_auc_folds > 1) auc[v] <- cv_auc(yA, X[u, , drop = FALSE], fit_fun, cv_auc_folds)
    minfo[[v]] <- data.frame(visit = v, n_active = length(a), n_unc = length(u), n_params = p,
                             n_A1 = sum(yA == 1), n_A0 = sum(yA == 0),
                             epv_A = min(sum(yA == 1), sum(yA == 0)) / p,
                             n_cens = sum(yC == 0),
                             epv_C = min(sum(yC == 1), sum(yC == 0)) / p)
  }
  model_info <- do.call(rbind, minfo)
  lowA <- model_info$visit[which(model_info$epv_A < epv_min)]
  lowC <- model_info$visit[which(model_info$epv_C < epv_min)]
  if (length(lowA))
    warning("Treatment model has < ", epv_min, " events (or non-events) per parameter at visits: ",
            paste(lowA, collapse = ","), ". Consider fit_fun = fit_glmnet or fewer covariates.")
  if (length(lowC))
    warning("Censoring model has < ", epv_min, " censored (or uncensored) women per parameter at visits: ",
            paste(lowC, collapse = ","), ". Consider fit_fun_C = fit_glmnet or fewer covariates.")
  
  # ---- cumulative weights for a rule and truncation bound ----
  compute_w <- function(ro, bound) {
    gA <- ifelse(ro$d == 1, pA, 1 - pA)            # P(A_v = d_v | history)
    f <- ifelse(act, pmax(gA, bound) * pmax(pC, bound), 1)
    cumw <- t(apply(f, 1, cumprod))
    list(gA = gA, H = ro$I / cumw)                  # H_v = I_v / prod_{s<=v} g_A g_C
  }
  
  # ---- 1. adherence table ----
  adherence <- do.call(rbind, lapply(names(rule_obj), function(r) {
    ro <- rule_obj[[r]]
    do.call(rbind, lapply(1:K, function(v) {
      Iprev <- if (v == 1) rep(TRUE, n) else ro$I[, v - 1]
      enter <- act[, v] & Iprev
      cs <- cstar[, v + 1] == 1
      miss <- Dm[, v + 1] %in% 1
      req <- ro$d[, v] == 1                                   # rule requires MHT at v
      newly <- req & (if (v == 1) rep(TRUE, n) else ro$d[, v - 1] == 0)  # first required visit
      n_sr <- sum(enter & unc[, v] & newly)
      n_sf <- sum(enter & unc[, v] & newly & ro$I[, v])
      data.frame(
        rule = r, visit = v,
        n_enter = sum(enter),
        n_censored_other = sum(enter & cs & !miss),
        n_censored_missing_mht = sum(enter & cs & miss),
        n_deviate = sum(enter & unc[, v] & (A[, v + 1] != ro$d[, v]), na.rm = TRUE),
        n_follow = sum(act[, v] & ro$I[, v]),
        n_required = sum(enter & req),
        n_required_follow = sum(act[, v] & ro$I[, v] & req),
        n_start_required = n_sr,
        n_start_follow = n_sf,
        pct_start_follow = if (n_sr > 0) n_sf / n_sr else NA_real_,
        n_following_incl_events = sum(ro$I[, v]),
        pct_baseline_following = sum(ro$I[, v]) / n)
    }))
  }))
  
  # ---- 2. positivity ----
  positivity <- do.call(rbind, lapply(names(rule_obj), function(r) {
    ro <- rule_obj[[r]]
    gA <- ifelse(ro$d == 1, pA, 1 - pA)
    do.call(rbind, lapply(1:K, function(v) {
      Iprev <- if (v == 1) rep(TRUE, n) else ro$I[, v - 1]
      enter <- act[, v] & Iprev
      x <- gA[enter, v]; z <- pC[enter, v]
      row <- list(rule = r, visit = v, n_enter = sum(enter),
                  gA_min = mfun(x, min), gA_p01 = qfun(x, .01), gA_p05 = qfun(x, .05),
                  gA_p50 = qfun(x, .5), gC_min = mfun(z, min), gC_p05 = qfun(z, .05),
                  gC_p50 = qfun(z, .5), auc_A_cv = auc[v])
      for (b in bounds) {
        row[[paste0("pct_gA_lt_", b)]] <- mfun(x, function(q) mean(q < b))
        row[[paste0("pct_gC_lt_", b)]] <- mfun(z, function(q) mean(q < b))
      }
      as.data.frame(row, check.names = FALSE)
    }))
  }))
  
  # ---- 3. cumulative weights ----
  weights <- do.call(rbind, lapply(bounds, function(b) {
    do.call(rbind, lapply(names(rule_obj), function(r) {
      ro <- rule_obj[[r]]; w <- compute_w(ro, b)
      do.call(rbind, lapply(1:K, function(v) {
        sel <- act[, v] & ro$I[, v]
        h <- w$H[sel, v]
        data.frame(
          rule = r, visit = v, bound = b, n = length(h),
          ess = ess_fun(h), p50 = qfun(h, .5), p90 = qfun(h, .9), p99 = qfun(h, .99),
          max_w = mfun(h, max),
          max_share = mfun(h, function(q) max(q) / sum(q)),
          mean_H_all = mean(w$H[, v]),   # should be ~1 for small bounds if g is adequate
          pct_trunc = mfun(which(sel), function(i)
            mean(w$gA[i, v] < b | pC[i, v] < b)))
      }))
    }))
  }))
  
  # ---- 4. pooled events (rule-independent, not tabulated by treatment) ----
  events <- data.frame(
    visit = 1:K,
    n_at_risk = sapply(1:K, function(v) sum(act[, v])),
    n_censored = sapply(1:K, function(v) sum(act[, v] & cstar[, v + 1] == 1)),
    n_events = sapply(1:K, function(v) sum(ev_time == v)))
  events$cum_events <- cumsum(events$n_events)
  
  res <- list(adherence = adherence, positivity = positivity, weights = weights,
              events = events, model_info = model_info, n = n, K = K, bounds = bounds,
              rules = names(rules), covs = covs, hot_codes = if (hot_rule) hot_codes else NULL,
              tau_stage = tau_S, tau_hot = tau_H)
  
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    for (nm in c("adherence", "positivity", "weights", "events", "model_info"))
      write.csv(res[[nm]], file.path(outdir, paste0(nm, ".csv")), row.names = FALSE)
  }
  res
}

# ---------- data-adaptive K and delta ---------------------------------------
# thr is a list you PRE-SPECIFY before looking at the tables, e.g.
#   thr <- list(N = 30, ESS = 30, p99 = 20, max_share = 0.10, pct_trunc = 0.05, V = 20)
# K = largest visit such that, for every rule in `rules` and EVERY visit up to it,
# all criteria hold. With both rule families in `res`, pass `rules` to evaluate one
# family at a time, e.g. rules = c("never", "hot_early").

select_K <- function(res, thr, bound, rules = NULL) {
  w <- res$weights[abs(res$weights$bound - bound) < 1e-12, ]
  if (!is.null(rules)) w <- w[w$rule %in% rules, ]
  w <- merge(w, res$events[, c("visit", "cum_events")], by = "visit")
  w$pass <- with(w, n >= thr$N & ess >= thr$ESS & p99 <= thr$p99 &
                   max_share <= thr$max_share & pct_trunc <= thr$pct_trunc &
                   cum_events >= thr$V)
  w$pass[is.na(w$pass)] <- FALSE
  ok <- tapply(w$pass, w$visit, all)
  K_sel <- 0
  for (v in sort(as.integer(names(ok)))) {
    if (isTRUE(ok[[as.character(v)]])) K_sel <- v else break
  }
  list(K = K_sel, detail = w[order(w$rule, w$visit), ])
}

# Largest delta whose feasible K is at least delta + margin.
# Stage rules: early_rule = "early", late_prefix = "late_d" (defaults).
# Symptom rules: early_rule = "hot_early", late_prefix = "hot_d".
pick_delta <- function(res, thr, bound, deltas, margin = 2,
                       early_rule = "early", late_prefix = "late_d") {
  out <- data.frame(delta = deltas, K = NA_integer_)
  for (i in seq_along(deltas))
    out$K[i] <- select_K(res, thr, bound,
                         rules = c("never", early_rule, paste0(late_prefix, deltas[i])))$K
  out$ok <- out$K >= out$delta + margin
  list(table = out, delta = if (any(out$ok)) max(out$delta[out$ok]) else NA)
}

# ---------- plots ------------------------------------------------------------

plot_diagnostics <- function(res, bound, file = "diagnostics.pdf") {
  pdf(file, width = 9, height = 4)
  cols <- seq_along(res$rules); names(cols) <- res$rules
  w <- res$weights[abs(res$weights$bound - bound) < 1e-12, ]
  op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1)); on.exit({par(op); dev.off()})
  plot_one <- function(df, col, ylab, log = "") {
    plot(NA, xlim = c(1, res$K), ylim = range(df[[col]], na.rm = TRUE),
         xlab = "Visit", ylab = ylab, log = log)
    for (r in res$rules) {
      d <- df[df$rule == r, ]; lines(d$visit, d[[col]], col = cols[[r]], lwd = 2)
    }
  }
  plot_one(res$adherence, "pct_baseline_following", "Fraction following rule")
  legend("topright", res$rules, col = cols, lwd = 2, cex = .7, bty = "n")
  plot_one(w, "ess", "Effective sample size", log = "y")
  plot_one(w, "p99", "99th pct of cumulative weight", log = "y")
  invisible(NULL)
}

# ---------- example usage ----------------------------------------------------
# dat <- read.csv("mht_final.csv", check.names = FALSE)
# res <- run_diagnostics(dat, K = 10,
#          status_trigger = c(2, 1),     # late peri or natural post (stage rules)
#          deltas = c(1, 2),
#          hot_rule = TRUE, hot_codes = 1:4,   # any hot flashes; try 3:4 as a sensitivity run
#          exclude_status0 = c(0, 7, 5), # optional: BSO, hysterectomy, pregnant at baseline
#          fit_fun_C = fit_glmnet,       # optional: penalized censoring model
#          bounds = c(0.01, 0.05), outdir = "diag_out")
# res$model_info                         # check events-per-parameter by visit
# thr <- list(N = 30, ESS = 30, p99 = 20, max_share = 0.10, pct_trunc = 0.05, V = 20)
# select_K(res, thr, 0.01, rules = c("never", "hot_early"))$K

dat <- final_df
res <- run_diagnostics(dat, K = 10,
         status_trigger = c(2, 1),     # late peri or natural post (stage rules)
         deltas = c(1, 2),
         hot_rule = TRUE, hot_codes = 1:4,   # any hot flashes; try 3:4 as a sensitivity run
         exclude_status0 = c(0, 7, 5), # optional: BSO, hysterectomy, pregnant at baseline
         fit_fun_C = fit_glmnet,       # optional: penalized censoring model
         bounds = c(0.01, 0.05), outdir = "diag_out")
res$model_info                         # check events-per-parameter by visit
thr <- list(N = 30, ESS = 30, p99 = 20, max_share = 0.10, pct_trunc = 0.05, V = 20)
select_K(res, thr, 0.01, rules = c("never", "hot_early"))$K
source(here("02-adherence-frames.R"))
adherence_frames(res$adherence, file = here("presentation", "adherence_filled.tex"))
adherence_frames(res$adherence, standalone = FALSE, file = here("presentation", "frames.tex"))  # frames only, to \input into your deck
