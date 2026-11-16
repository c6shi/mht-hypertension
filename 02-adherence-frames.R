# 02-adherence-frames.R
# Writes one Beamer frame per rule, showing that rule's adherence table.
#   adherence_frames(file = "adherence_template.tex")                 # blank cells to fill by hand
#   adherence_frames(res$adherence, file = "adherence_filled.tex")    # filled from run_diagnostics()
#   adherence_frames(res$adherence, standalone = FALSE, file = "frames.tex")  # frames only, for \input
# Columns shown (from the adherence table):
#   Entering   = n_enter                  event-free, uncensored and following the rule through the previous visit
#   Other      = n_censored_other         censored for a reason other than missing MHT
#   Missing    = n_censored_missing_mht   censored because MHT was missing at the visit
#   Deviated   = n_deviate                MHT use differed from what the rule requires
#   Followed   = n_follow                 followed the rule (includes women whose event occurs at that visit)
#   Must start = n_start_required         uncensored women the rule requires to START MHT at this visit
#   Started    = n_start_follow (pct_start_follow)
# Identity (holds in the data): Entering = Other + Missing + Deviated + Followed.

adherence_frames <- function(adh = NULL,
                             rules = if (is.null(adh)) c("never", "early", "late_d1", "late_d2", "hot_early")
                             else unique(adh$rule),
                             file = "adherence_slides.tex", K = 10, standalone = TRUE) {
  titles <- c(never = "Never initiate",
              early = "Early initiation (stage-triggered)",
              late_d1 = "Delayed initiation: 1 visit after stage trigger",
              late_d2 = "Delayed initiation: 2 visits after stage trigger",
              hot_early = "Early initiation (symptom-triggered)")
  descr <- c(never = "MHT is never started; any MHT use is a deviation.",
             early = "Start MHT in the interval after first reaching late perimenopause or later, then continue.",
             late_d1 = "Start MHT one visit after first reaching late perimenopause or later, then continue.",
             late_d2 = "Start MHT two visits after first reaching late perimenopause or later, then continue.",
             hot_early = "Start MHT in the interval after first reporting hot flashes, then continue.")
  esc <- function(x) gsub("_", "\\\\_", x)
  
  int <- function(x) ifelse(is.na(x), "--", format(as.integer(x), trim = TRUE))
  frame_for <- function(r) {
    ttl <- if (r %in% names(titles)) titles[[r]] else esc(r)
    dsc <- if (r %in% names(descr)) descr[[r]] else ""
    rows <- character(K)
    if (is.null(adh)) {
      for (v in 1:K) rows[v] <- sprintf("%d & & & & & & & \\\\", v)
    } else {
      d <- adh[adh$rule == r, ]
      if (nrow(d) == 0) stop("rule not found in adherence table: ", r)
      d <- d[order(d$visit), ]
      for (v in 1:K) {
        x <- d[d$visit == v, ]
        if (nrow(x) == 0) { rows[v] <- sprintf("%d & -- & -- & -- & -- & -- & -- & -- \\\\", v); next }
        started <- if (is.na(x$n_start_required) || x$n_start_required == 0) "--" else
          sprintf("%d (%.0f\\%%)", as.integer(x$n_start_follow), 100 * x$pct_start_follow)
        rows[v] <- paste0(v, " & ", int(x$n_enter), " & ", int(x$n_censored_other), " & ",
                          int(x$n_censored_missing_mht), " & ", int(x$n_deviate), " & ",
                          int(x$n_follow), " & ", int(x$n_start_required), " & ", started, " \\\\")
      }
    }
    c(sprintf("\\begin{frame}{Adherence: %s}", ttl),
      "\\centering",
      sprintf("{\\small \\textbf{Rule:} %s}\\par\\vspace{6pt}", dsc),
      "{\\scriptsize",
      "\\begin{tabular}{r r r r r r r r}",
      "\\toprule",
      " & & \\multicolumn{2}{c}{Censored} & & & & \\\\",
      "\\cmidrule(lr){3-4}",
      "Visit & Entering & Other & Missing MHT & Deviated & Followed & Must start & Started (\\%) \\\\",
      "\\midrule", rows, "\\bottomrule", "\\end{tabular}}",
      "\\par\\vspace{6pt}",
      "{\\tiny Entering = event-free, uncensored and following the rule through the previous visit; Entering = Censored + Deviated + Followed. Started = women required to start at this visit who did.}",
      "\\end{frame}", "")
  }
  
  body <- unlist(lapply(rules, frame_for))
  out <- if (standalone)
    c("\\documentclass[aspectratio=169]{beamer}", "\\usepackage{booktabs}",
      "\\setbeamertemplate{navigation symbols}{}", "\\begin{document}", "", body, "\\end{document}")
  else body
  writeLines(out, file)
  invisible(file)
}
