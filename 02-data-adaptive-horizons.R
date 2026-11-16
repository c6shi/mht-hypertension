# Data adaptive horizon: choosing the best K (outcome blind)!
# However, we can choose to report all K's if we have time.

library(here)
library(dplyr)
library(tidyverse)
library(data.table)

final_df <- read.csv(here("data", "mht_final.csv"))

# Metrics for each horizon K: 
# 1) Adherence table
# 2) Positivity check
# 3) Cumulative-weight diagnostics (clever covariate)

# WRITTEN BY CLAUDE

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
# Aliased columns (e.g., empty factor levels) get coefficient 0.
fit_glm <- function(y, X_fit, X_new) {
  if (length(unique(y)) < 2) return(rep(mean(y), nrow(X_new)))
  f <- suppressWarnings(glm.fit(X_fit, y, family = binomial()))
  b <- f$coefficients; b[is.na(b)] <- 0
  as.numeric(plogis(X_new %*% b))
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

# ---------- main function ----------------------------------------------------

run_diagnostics <- function(dat, K = 10, status_trigger, deltas = c(1, 2),
                            bounds = c(0.01, 0.05),
                            rhs = "AGE + STATUS_prev + HOT_prev + MHT_prev + ever_trig_prev",
                            fit_fun = fit_glm, cv_auc_folds = 5,
                            outdir = NULL, seed = 1) {
  if (missing(status_trigger))
    stop("Set status_trigger to the STATUS codes meaning late peri-/early post-menopause.")
  set.seed(seed)
  
  # new-user, hypertension-free at baseline
  dat <- dat[which(dat$MHT0 %in% 0 & dat$HTN0 %in% 0), , drop = FALSE]
  n <- nrow(dat)
  if (anyNA(dat$AGE)) stop("AGE has missing values: impute or drop before running.")
  message("Analytic sample (MHT-naive, HTN-free at baseline): n = ", n)
  
  S  <- get_mat(dat, "STATUS", K);  Hf <- get_mat(dat, "HOTFLAS", K)
  A  <- get_mat(dat, "MHT", K);     Dm <- get_mat(dat, "D_MHT", K)
  Cc <- get_mat(dat, "C_HTN", K);   Y  <- get_mat(dat, "HTN", K)
  A[, 1] <- 0
  
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
  
  # trigger visit
  trig <- matrix(FALSE, n, K + 1)
  for (v in 0:K) trig[, v + 1] <- S[, v + 1] %in% status_trigger
  lim <- pmin(ev_time, cens_time - 1, K)
  tau <- rep(Inf, n)
  for (v in rev(0:(K - 1))) tau[which(trig[, v + 1] & (v <= lim))] <- v
  
  # ---- rules ----
  rules <- list(never = NA, early = 0)
  for (dl in deltas) rules[[paste0("late_d", dl)]] <- dl
  
  make_d <- function(delta) {
    M <- matrix(0, n, K)
    if (!is.na(delta)) for (v in 1:K) M[, v] <- as.numeric(tau <= v - 1 - delta)
    M
  }
  
  rule_obj <- lapply(rules, function(delta) {
    d <- make_d(delta)
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
  
  # ---- treatment and censoring mechanisms (rule-independent) ----
  lev_S <- sort(unique(S[!is.na(S)])); lev_H <- sort(unique(Hf[!is.na(Hf)]))
  fac <- function(x, lev)
    factor(ifelse(is.na(x), "NA", as.character(x)), levels = c(as.character(lev), "NA"))
  design <- function(v) {
    Xdf <- data.frame(
      AGE = dat$AGE,
      STATUS_prev = fac(S[, v], lev_S),     # column v = visit v-1
      HOT_prev = fac(Hf[, v], lev_H),
      MHT_prev = ifelse(is.na(A[, v]), 0, A[, v]),
      ever_trig_prev = as.numeric(tau <= v - 1))
    model.matrix(as.formula(paste("~", rhs)), data = Xdf)
  }
  
  pA <- matrix(1, n, K); pC <- matrix(1, n, K); auc <- rep(NA_real_, K)
  for (v in 1:K) {
    a <- which(act[, v]); u <- a[unc[a, v]]
    if (length(a) < 20 || length(u) < 20) {
      warning("Visit ", v, ": fewer than 20 women at risk; g set to 1. Interpret nothing past here.")
      next
    }
    X <- design(v)
    pC[a, v] <- fit_fun(1 - cstar[a, v + 1], X[a, , drop = FALSE], X[a, , drop = FALSE])
    yA <- A[u, v + 1]
    pA[a, v] <- fit_fun(yA, X[u, , drop = FALSE], X[a, , drop = FALSE])
    if (cv_auc_folds > 1) auc[v] <- cv_auc(yA, X[u, , drop = FALSE], fit_fun, cv_auc_folds)
  }
  
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
      data.frame(
        rule = r, visit = v,
        n_enter = sum(enter),
        n_censored_other = sum(enter & cs & !miss),
        n_censored_missing_mht = sum(enter & cs & miss),
        n_deviate = sum(enter & unc[, v] & (A[, v + 1] != ro$d[, v]), na.rm = TRUE),
        n_follow = sum(act[, v] & ro$I[, v]),
        n_triggered_enter = sum(enter & tau <= v - 1),
        n_triggered_follow = sum(act[, v] & ro$I[, v] & tau <= v - 1),
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
              events = events, n = n, K = K, bounds = bounds,
              rules = names(rules), tau = tau)
  
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    for (nm in c("adherence", "positivity", "weights", "events"))
      write.csv(res[[nm]], file.path(outdir, paste0(nm, ".csv")), row.names = FALSE)
  }
  res
}

# ---------- data-adaptive K and delta ---------------------------------------
# thr is a list you PRE-SPECIFY before looking at the tables, e.g.
#   thr <- list(N = 30, ESS = 30, p99 = 20, max_share = 0.10, pct_trunc = 0.05, V = 20)
# K = largest visit such that, for every rule and EVERY visit up to it, all criteria hold.

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
pick_delta <- function(res, thr, bound, deltas, margin = 2) {
  out <- data.frame(delta = deltas, K = NA_integer_)
  for (i in seq_along(deltas))
    out$K[i] <- select_K(res, thr, bound,
                         rules = c("never", "early", paste0("late_d", deltas[i])))$K
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

# ---------- example: SuperLearner learner -----------------------------------
# fit_sl <- function(y, X_fit, X_new) {
#   library(SuperLearner)
#   Xf <- as.data.frame(X_fit[, -1, drop = FALSE])   # drop intercept column
#   Xn <- as.data.frame(X_new[, -1, drop = FALSE])
#   if (length(unique(y)) < 2) return(rep(mean(y), nrow(Xn)))
#   sl <- SuperLearner(Y = y, X = Xf, newX = Xn, family = binomial(),
#                      SL.library = c("SL.glm", "SL.mean", "SL.gam"))
#   as.numeric(sl$SL.predict)
# }
# res <- run_diagnostics(dat, K = 10, status_trigger = c(...), fit_fun = fit_sl)

dat <- final_df
res <- run_diagnostics(dat, K = 10, status_trigger = c(3, 2), deltas = c(1, 2, 3),
                       bounds = c(0.001, 0.01, 0.05), outdir = "diag_out")
options(width = 200)
cat("\n--- adherence (early, late_d2) ---\n")
print(subset(res$adherence, rule %in% c("early", "late_d2")))
cat("\n--- identity n_enter = cens + deviate + follow ---\n")
a <- res$adherence
print(all(a$n_enter == a$n_censored_other + a$n_censored_missing_mht + a$n_deviate + a$n_follow))
cat("\n--- weights (bound 0.01) ---\n")
print(subset(res$weights, bound == 0.01 & rule %in% c("early", "late_d2")), digits = 3)
cat("\n--- mean_H_all by rule, bound 0.001 (should be near 1) ---\n")
w <- res$weights[res$weights$bound == 0.001, ]
print(round(tapply(w$mean_H_all, list(w$visit, w$rule), mean), 2))
cat("\n--- events ---\n"); print(res$events)
thr <- list(N = 30, ESS = 30, p99 = 20, max_share = 0.10, pct_trunc = 0.05, V = 20)
k <- select_K(res, thr, bound = 0.01); cat("\nK selected:", k$K, "\n")
print(pick_delta(res, thr, 0.01, c(1, 2, 3)))
plot_diagnostics(res, 0.01, "diag_out/diagnostics.pdf")
