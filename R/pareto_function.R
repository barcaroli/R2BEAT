# Pareto frontier allocation: adjusts a global scaling factor t on domain/variable CV
# targets, computes allocations via R2BEAT, and finds the smallest feasible total
# sample size within the target/tolerance band using bracketing + bisection, while
# honoring optional CV caps. The resulting solution is the closest-to-target
# feasible allocation on the Pareto frontier between precision and cost.
#
# Parameters:
#   strata         Data frame with strata information. Required columns: STRATUM,
#                  N, COST, CENS (case-insensitive).
#   current_cvs    Data frame with initial CV targets per domain. Must include DOM
#                  plus CV1..CVp (or V1..Vp, which are renamed to CV*).
#   target_size    Desired total sample size (scalar).
#   tolerance      Allowed shortfall from target_size (scalar, >= 0).
#   cv_caps        Optional data frame with columns DOM, VAR, MAX_CV to cap CVs
#                  by domain/variable. VAR can be V#, CV#, or integer index.
#   t0             Initial scaling factor for CV targets.
#   inactive_cv    CV values >= inactive_cv are treated as inactive targets.
#   min_cv         Lower bound for active CV targets after scaling.
#   eps_inact      Small tolerance when detecting inactive CVs.
#   bracket_fac    Multiplicative step used to expand the bracketing interval.
#   max_iter       Maximum bisection iterations.
#   max_same_iter  Stop bracketing if n(t) repeats this many times (caps plateau).
#   minnumstrat    Minimum sample size per stratum passed to beat1cv_fun.
#   beat1cv_fun    Optional function to compute expected CVs (e.g., beat.1cv_new).
#   verbose        If TRUE, print progress messages.
#   plot           If TRUE, plot CVs by domain and allocations.
#   plot_dir       Output directory for plots.
#   plot_prefix    Filename prefix for plot outputs.
#   save_pdf       If TRUE, save plots to a PDF file.
#   show_targets   If TRUE, overlay initial target CVs on plots.
#   show_caps      If TRUE, overlay CV caps on plots (when provided).
#   plot_convergence If TRUE, create convergence plot showing bracketing and bisection
pareto <- function(strata,
                   current_cvs,
                   target_size,
                   tolerance,
                   cv_caps = NULL,         
                   t0 = 1,
                   inactive_cv = 0.999,
                   min_cv = 1e-8,
                   eps_inact = 1e-12,
                   bracket_fac = 1.25,
                   max_iter = 50,
                   max_same_iter = max_iter,  
                   minnumstrat = 2,
                   beat1cv_fun = NULL,      
                   verbose = TRUE,
                   plot = FALSE,
                   plot_dir = "outputs",
                   plot_prefix = "pareto_run",
                   save_pdf = TRUE,
                   show_targets = TRUE,     
                   show_caps = FALSE,
                   plot_convergence = TRUE  
) {
  
  # ---- checks ----
  if (!requireNamespace("R2BEAT", quietly = TRUE)) stop("Package R2BEAT is required.")
  if (!is.data.frame(strata)) stop("'strata' must be a dataframe.")
  if (!is.data.frame(current_cvs)) stop("'current_cvs' must be a dataframe.")
  if (!is.numeric(target_size) || length(target_size) != 1) stop("'target_size' must be a numeric scalar.")
  if (!is.numeric(tolerance)  || length(tolerance)  != 1) stop("'tolerance' must be a numeric scalar.")
  if (tolerance < 0) stop("'tolerance' deve essere >= 0.")
  
  if (!is.numeric(max_same_iter) || length(max_same_iter) != 1 || max_same_iter < 1) {
    stop("'max_same_iter' deve essere un intero >= 1.")
  }
  names(strata) <- toupper(names(strata))
  names(current_cvs) <- toupper(names(current_cvs))
  
  if ("C" %in% names(strata)) names(strata)[names(strata) == "C"] <- "COST"
  if ("STRATO" %in% names(strata)) names(strata)[names(strata) == "STRATO"] <- "STRATUM"
  
  req_cols <- c("STRATUM", "N", "COST", "CENS")
  miss <- setdiff(req_cols, names(strata))
  if (length(miss) > 0) stop("In 'strata' the following columns are missing: ", paste(miss, collapse=", "))
  
  if (!("DOM" %in% names(current_cvs))) stop("In 'current_cvs' column DOM is missing.")
  
  # Accept CV1..CVp or V1..Vp
  if (length(grep("^CV[0-9]+$", names(current_cvs))) == 0) {
    vcols <- grep("^V[0-9]+$", names(current_cvs), value = TRUE)
    if (length(vcols) == 0) stop("I can't find neither CV1..CVp nor V1..Vp in 'current_cvs'.")
    names(current_cvs)[match(vcols, names(current_cvs))] <- sub("^V", "CV", vcols)
  }
  
  cv_cols <- grep("^CV[0-9]+$", names(current_cvs), value = TRUE)
  cv_cols <- cv_cols[order(as.integer(sub("^CV", "", cv_cols)))]
  n_var <- length(cv_cols)
  vcols_expected <- paste0("V", 1:n_var)
  
  # Initial targets (for plot)
  errors_initial <- current_cvs[, c("DOM", cv_cols), drop = FALSE]
  
  vlog <- function(...) if (isTRUE(verbose)) cat(sprintf(...), "\n")
  in_band <- function(n) (n <= target_size) && (n >= (target_size - tolerance))
  
  # ------------------------------------------------------------
  # NEW: Initialize tracking for convergence plot
  # ------------------------------------------------------------
  convergence_history <- list(
    t_vals = numeric(0),
    n_vals = numeric(0),
    phase = character(0),
    iteration = integer(0)
  )
  
  # ------------------------------------------------------------
  # Caps: build a caps matrix (DOM x CVk) with NA where missing
  # ------------------------------------------------------------
  caps_mat <- NULL
  if (!is.null(cv_caps)) {
    if (!all(c("DOM","VAR","MAX_CV") %in% toupper(names(cv_caps)))) {
      stop("cv_caps must have columns: DOM, VAR, MAX_CV (case-insensitive).")
    }
    cv_caps2 <- cv_caps
    names(cv_caps2) <- toupper(names(cv_caps2))
    cv_caps2$DOM <- as.character(cv_caps2$DOM)
    
    # Normalize VAR -> index 1..n_var
    var_to_idx <- function(v) {
      if (is.numeric(v)) return(as.integer(v))
      v <- toupper(as.character(v))
      v <- sub("^CV", "", v)
      v <- sub("^V",  "", v)
      as.integer(v)
    }
    cv_caps2$IDX <- var_to_idx(cv_caps2$VAR)
    if (any(is.na(cv_caps2$IDX)) || any(cv_caps2$IDX < 1) || any(cv_caps2$IDX > n_var)) {
      stop("cv_caps$VAR must identify a variable among 1..", n_var, " (ex. V3, CV3 or 3).")
    }
    
    # Caps matrix aligned to domains in errors_initial
    dom_levels <- as.character(errors_initial$DOM)
    caps_mat <- matrix(NA_real_, nrow = length(dom_levels), ncol = n_var,
                       dimnames = list(dom_levels, paste0("CV", 1:n_var)))
    
    for (i in seq_len(nrow(cv_caps2))) {
      d <- cv_caps2$DOM[i]
      k <- cv_caps2$IDX[i]
      m <- cv_caps2$MAX_CV[i]
      if (!is.finite(m) || m <= 0 || m >= inactive_cv) {
        stop("MAX_CV must be in (0, ", inactive_cv, "). Found: ", m)
      }
      if (!(d %in% dom_levels)) {
        stop("DOM in cv_caps not found in current_cvs: ", d)
      }
      caps_mat[d, paste0("CV", k)] <- m
    }
  }
  
  # ------------------------------------------------------------
  # Build "base targets" used by optimization:
  #   base = min(initial_target, cap) when cap is present
  # ------------------------------------------------------------
  errors_base <- errors_initial
  if (!is.null(caps_mat)) {
    for (k in 1:n_var) {
      ck <- paste0("CV", k)
      cap_col <- caps_mat[errors_base$DOM, ck]
      x <- errors_base[[ck]]
      has_cap <- !is.na(cap_col)
      # If cap is present, enforce target <= cap (even if previously 0.999)
      x[has_cap] <- pmin(x[has_cap], cap_col[has_cap])
      errors_base[[ck]] <- x
    }
  }
  
  # ------------------------------------------------------------
  # Scaling for the dual with caps:
  #   scaled = base * t for active cells
  #   then clamp to [min_cv, inactive_cv] and also clamp to cap (if present)
  # ------------------------------------------------------------
  scale_errors <- function(base_df, t) {
    out <- base_df
    for (k in 1:n_var) {
      ck <- paste0("CV", k)
      x <- out[[ck]]
      
      active <- !is.na(x) & (x < (inactive_cv - eps_inact))
      x_new <- x
      x_new[active] <- pmin(inactive_cv, pmax(min_cv, x[active] * t))
      
      # Clamp to cap (maximum allowed) if present
      if (!is.null(caps_mat)) {
        cap_col <- caps_mat[out$DOM, ck]
        has_cap <- !is.na(cap_col)
        x_new[has_cap] <- pmin(x_new[has_cap], cap_col[has_cap])
      }
      
      out[[ck]] <- x_new
    }
    out
  }
  
  alloc_for_t <- function(t) {
    e_t <- scale_errors(errors_base, t)
    a <- R2BEAT::beat.1st(stratif = strata, errors = e_t)
    # beat.1st alloc table includes a total row; drop it to match strata rows
    n_vec <- a$alloc$ALLOC
    if (length(n_vec) == (nrow(strata) + 1L)) n_vec <- n_vec[-length(n_vec)]
    list(errors = e_t, alloc = a, n_tot = sum(n_vec))
  }
  n_for_t <- function(t) alloc_for_t(t)$n_tot
  
  # ------------------------------------------------------------
  # NEW: Helper function to track convergence
  # ------------------------------------------------------------
  track_point <- function(t, n, phase, iter) {
    convergence_history$t_vals <<- c(convergence_history$t_vals, t)
    convergence_history$n_vals <<- c(convergence_history$n_vals, n)
    convergence_history$phase <<- c(convergence_history$phase, phase)
    convergence_history$iteration <<- c(convergence_history$iteration, iter)
  }
  
  # ------------------------------------------------------------
  # Bracket + bisection (with possible plateau due to caps)
  # ------------------------------------------------------------
  n0 <- n_for_t(t0)
  vlog("Start: t=%.10f -> n=%d", t0, n0)
  track_point(t0, n0, "Start", 0)  # NEW
  
  # If n0 > target_size: need to relax (increase t)
  # If n0 <= target_size: we can tighten (decrease t)
  if (n0 > target_size) {
    t_lo <- t0; n_lo <- n0
    t_hi <- t0; n_hi <- n0
    same_n_hi <- 0L
    prev_n_hi <- NA_integer_
    bracket_iter <- 0  # NEW
    repeat {
      t_hi <- t_hi * bracket_fac
      n_hi <- n_for_t(t_hi)
      bracket_iter <- bracket_iter + 1  # NEW
      vlog("Bracket-up: t=%.10f -> n=%d", t_hi, n_hi)
      track_point(t_hi, n_hi, "Bracket-up", bracket_iter)  # NEW
      
      # Check if n_hi has plateaued
      if (!is.na(prev_n_hi) && (n_hi == prev_n_hi)) {
        same_n_hi <- same_n_hi + 1L
        if (same_n_hi >= max_same_iter) {
          stop("Bracketing failed: n(t) is stuck at ", n_hi, " for ", max_same_iter,
               " iterations (infeasible due to caps or other constraints).")
        }
      } else {
        same_n_hi <- 0L
      }
      prev_n_hi <- n_hi
      
      if (n_hi <= target_size) {
        vlog("  --> Bracket found (t_lo=%.10f infeasible, t_hi=%.10f feasible)", t_lo, t_hi)
        break
      }
    }
  } else {
    t_hi <- t0; n_hi <- n0
    t_lo <- t0; n_lo <- n0
    same_n_lo <- 0L
    prev_n_lo <- NA_integer_
    bracket_iter <- 0  # NEW
    repeat {
      t_lo <- t_lo / bracket_fac
      n_lo <- n_for_t(t_lo)
      bracket_iter <- bracket_iter + 1  # NEW
      vlog("Bracket-down: t=%.10f -> n=%d", t_lo, n_lo)
      track_point(t_lo, n_lo, "Bracket-down", bracket_iter)  # NEW
      
      # Check plateau
      if (!is.na(prev_n_lo) && (n_lo == prev_n_lo)) {
        same_n_lo <- same_n_lo + 1L
        if (same_n_lo >= max_same_iter) {
          res <- alloc_for_t(t_hi)
          alloc_df <- data.frame(STRATUM = strata$STRATUM, n = res$alloc$n)
          
          # Return infeasible result with tracking
          return(structure(list(
            feasible = FALSE,
            message = sprintf("Can't relax below n=%d (same for %d iterations in bracketing).", n_lo, max_same_iter),
            t = t_hi,
            n_total = n_hi,
            errors_scaled = res$errors,
            errors_initial = errors_initial,
            errors_base = errors_base,
            allocation = alloc_df,
            expectedCV = NULL,
            expectedCV_spec = NULL,
            binding = NULL,
            plot_pdf = NULL,
            convergence_history = convergence_history  # NEW
          )))
        }
      } else {
        same_n_lo <- 0L
      }
      prev_n_lo <- n_lo
      
      if (n_lo > target_size) {
        vlog("  --> Bracket found (t_lo=%.10f infeasible, t_hi=%.10f feasible)", t_lo, t_hi)
        break
      }
    }
  }
  
  # If already in tolerance band after bracketing
  if (in_band(n_hi)) {
    res <- alloc_for_t(t_hi)
    n_vec <- res$alloc$alloc$ALLOC
    if (length(n_vec) == (nrow(strata) + 1L)) n_vec <- n_vec[-length(n_vec)]
    alloc_df <- data.frame(STRATUM = strata$STRATUM, n = n_vec)
    
    # ---- expected CV ----
    expectedCV <- NULL
    expectedCV_spec <- NULL
    
    if (is.null(beat1cv_fun) && exists("beat.1cv_2", mode="function")) beat1cv_fun <- get("beat.1cv_2")
    if (!is.null(beat1cv_fun)) {
      n_vec <- res$alloc$alloc$ALLOC
      if (length(n_vec) == (nrow(strata) + 1L)) n_vec <- n_vec[-length(n_vec)]
      cv_out <- beat1cv_fun(stratif = strata, alloc = n_vec, minnumstrat = minnumstrat)
      expectedCV <- cv_out$expectedCV
      expectedCV_spec <- cv_out$expectedCV_spec
    }
    
    return(structure(list(
      feasible = TRUE,
      t = t_hi,
      n_total = res$n_tot,
      early_stop = TRUE,
      errors_scaled = res$errors,
      errors_initial = errors_initial,
      errors_base = errors_base,
      allocation = alloc_df,
      expectedCV = expectedCV,
      expectedCV_spec = expectedCV_spec,
      binding = NULL,
      plot_pdf = NULL,
      convergence_history = convergence_history  # NEW
    )))
  }
  
  # Now we have (t_lo infeasible) and (t_hi feasible)
  if (!(n_lo > target_size && n_hi <= target_size)) {
    stop("Bracket not valid (with caps): expected n_lo > target_size and n_hi <= target_size.")
  }
  
  best_t <- t_hi
  best_n <- n_hi
  early_stop <- in_band(best_n)
  
  if (!early_stop) {
    for (iter in 1:max_iter) {
      t_mid <- 0.5 * (t_lo + t_hi)
      n_mid <- n_for_t(t_mid)
      vlog("Bisect %02d: t=%.12f -> n=%d", iter, t_mid, n_mid)
      track_point(t_mid, n_mid, "Bisect", iter)  # NEW
      
      if (n_mid <= target_size) {
        t_hi <- t_mid
        if (n_mid >= best_n) {
          best_t <- t_mid
          best_n <- n_mid
        }
        if (in_band(n_mid)) { best_t <- t_mid; best_n <- n_mid; early_stop <- TRUE; break }
      } else {
        t_lo <- t_mid
      }
      
      if (abs(t_hi - t_lo) / max(1e-12, t_hi) < 1e-8) break
    }
  }
  
  res <- alloc_for_t(best_t)
  n_vec <- res$alloc$alloc$ALLOC
  if (length(n_vec) == (nrow(strata) + 1L)) n_vec <- n_vec[-length(n_vec)]
  alloc_df <- data.frame(STRATUM = strata$STRATUM, n = n_vec)
  
  # ---- expected CV ----
  expectedCV <- NULL
  expectedCV_spec <- NULL
  binding_tbl <- NULL
  
  if (is.null(beat1cv_fun) && exists("beat.1cv_2", mode="function")) beat1cv_fun <- get("beat.1cv_2")
  if (!is.null(beat1cv_fun)) {
    n_vec <- res$alloc$alloc$ALLOC
    if (length(n_vec) == (nrow(strata) + 1L)) n_vec <- n_vec[-length(n_vec)]
    cv_out <- beat1cv_fun(stratif = strata, alloc = n_vec, minnumstrat = minnumstrat)
    expectedCV <- cv_out$expectedCV
    expectedCV_spec <- cv_out$expectedCV_spec
  }
  
  # ------------------------------------------------------------
  # PLOTS (optional): CV by domain + allocations + CONVERGENCE
  #  - dashed line = INITIAL TARGETS (as requested)
  #  - optional second line = caps
  # ------------------------------------------------------------
  plot_pdf_path <- NULL
  
  plot_cv_domain <- function(dom_label, cv_vals, tar_vals = NULL, cap_vals = NULL) {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = c(8, 4, 4, 2) + 0.1)
    
    ymax <- max(cv_vals,
                if (!is.null(tar_vals)) tar_vals else 0,
                if (!is.null(cap_vals)) cap_vals else 0,
                na.rm = TRUE) * 1.15
    
    bp <- barplot(cv_vals,
                  names.arg = names(cv_vals),
                  las = 2, cex.names = 0.8,
                  col = "orange",
                  ylim = c(0, ymax),
                  main = paste("Expected CV by variable - Domain:", dom_label),
                  ylab = "CV")
    
    leg <- c("Achieved CV"); pch <- c(15); lty <- c(NA)
    
    if (!is.null(tar_vals)) {
      lines(bp, tar_vals, type="b", lty=2, pch=16)
      leg <- c(leg, "Initial target CV"); pch <- c(pch, 16); lty <- c(lty, 2)
    }
    if (!is.null(cap_vals)) {
      lines(bp, cap_vals, type="b", lty=3, pch=17)
      leg <- c(leg, "CV caps (max)"); pch <- c(pch, 17); lty <- c(lty, 3)
    }
    
    legend("topright", legend = leg, pch = pch, lty = lty, bty = "n", cex=0.6)
  }
  
  plot_alloc <- function(alloc_df_local) {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = c(9, 5, 4, 2) + 0.1)
    ymax <- max(alloc_df_local$n, na.rm = TRUE) * 1.15
    barplot(alloc_df_local$n,
            names.arg = alloc_df_local$STRATUM,
            las = 2, cex.names = 0.7,
            col = "orange",
            ylim = c(0, ymax),
            ylab = "Allocated sample size (n)",
            main = sprintf("Allocation by stratum (n_total=%d, target=%d)", sum(alloc_df_local$n), target_size))
  }
  
  # ------------------------------------------------------------
  # NEW: Convergence plot function
  # ------------------------------------------------------------
  plot_convergence_trajectory <- function(history, target, tol) {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = c(5, 5, 4, 2) + 0.1)
    
    # Determine plot limits
    xlim <- range(history$t_vals, na.rm = TRUE)
    ylim <- range(c(history$n_vals, target, target - tol), na.rm = TRUE)
    ylim <- ylim + c(-0.05, 0.05) * diff(ylim)
    xlim <- xlim + c(-0.04, 0.04) * diff(xlim)
    
    x_pad <- 0.02 * diff(xlim)
    x_left <- xlim[1] + x_pad
    x_right <- xlim[2] - x_pad
    clamp_center_x <- function(x, label, cex = 1) {
      half <- strwidth(label, cex = cex) / 2
      max_x <- x_right - half
      min_x <- x_left + half
      min(max(x, min_x), max_x)
    }
    clamp_right_x <- function(x, label, cex = 1) {
      max_x <- x_right - strwidth(label, cex = cex)
      min(max(x, x_left), max_x)
    }
    
    # Create base plot
    plot(history$t_vals, history$n_vals, type="n",
         xlim = xlim, ylim = ylim,
         xlab = "Scaling Factor (t)",
         ylab = "Sample Size (n)",
         main = "Pareto Bisection Search for Feasibility",
         cex.lab = 1.1, cex.main = 1.2)
    
    # Add budget line
    abline(h = target, col = "steelblue", lwd = 2, lty = 1)
    budget_x <- clamp_center_x(xlim[1] + 0.95 * diff(xlim), " ", cex = 0.7)
    text(budget_x, target, " ", pos = 3, col = "steelblue", cex = 0.7)
    
    # Add tolerance band (if any)
    if (tol > 0) {
      rect(xlim[1], target - tol, xlim[2], target, 
           col = rgb(0.7, 0.9, 1, 0.2), border = NA)
    }
    
    # Separate points by phase
    phases <- unique(history$phase)
    colors <- c("Start" = "red", "Bracket-up" = "pink", "Bracket-down" = "pink", 
                "Bisect" = "gray50")
    
    # Plot trajectory lines
    for (phase in phases) {
      idx <- which(history$phase == phase)
      if (length(idx) > 0) {
        points(history$t_vals[idx], history$n_vals[idx], 
               col = colors[phase], pch = 16, cex = 1.2)
        if (length(idx) > 1) {
          lines(history$t_vals[idx], history$n_vals[idx], 
                col = colors[phase], lwd = 1.5)
        }
      }
    }
    
    # Highlight key points
    # Start point
    start_idx <- which(history$phase == "Start")
    if (length(start_idx) > 0) {
      points(history$t_vals[start_idx], history$n_vals[start_idx], 
             pch = 21, cex = 2, bg = "red", col = "darkred", lwd = 2)
      start_label <- sprintf("Start (Infeasible)")
      start_x <- clamp_center_x(history$t_vals[start_idx], start_label, cex = 0.75)
      text(start_x, history$n_vals[start_idx], 
           start_label,
           pos = 3, cex = 0.75, col = "red", font = 2)
    }
    
    # Bracket found point (first feasible)
    bracket_idx <- which(history$phase %in% c("Bracket-up", "Bracket-down") & 
                           history$n_vals <= target)
    if (length(bracket_idx) > 0) {
      first_feas <- bracket_idx[1]
      points(history$t_vals[first_feas], history$n_vals[first_feas], 
             pch = 21, cex = 2, bg = "lightblue", col = "steelblue", lwd = 2)
      # bracket_label <- sprintf("Bracket found\n(Feasible)")
      bracket_label <- sprintf("Bracket found")
      bracket_x <- clamp_center_x(history$t_vals[first_feas], bracket_label, cex = 0.75)
      text(bracket_x, history$n_vals[first_feas], 
           bracket_label,
           pos = 1, cex = 0.75, col = "steelblue", font = 2)
    }
    
    # Converged point (final)
    final_idx <- length(history$t_vals)
    points(history$t_vals[final_idx], history$n_vals[final_idx], 
           pch = 21, cex = 2.5, bg = "darkblue", col = "black", lwd = 2)
    # converged_label <- sprintf("Converged: t=%.5f (n=%d)", history$t_vals[final_idx], history$n_vals[final_idx])
    converged_label <- sprintf("Converged                   ", history$t_vals[final_idx], history$n_vals[final_idx])
    converged_x <- clamp_right_x(history$t_vals[final_idx], converged_label, cex = 0.75)
    text(converged_x, history$n_vals[final_idx], 
         converged_label,
         pos = 4, cex = 0.75, col = "darkblue", font = 2)
    
    # Add annotations for infeasible/feasible regions
    mid_y <- mean(ylim)
    if (n0 > target) {
      # Infeasible region is at lower t values
      text(xlim[1] + 0.05*diff(xlim), ylim[2] - 0.1*diff(ylim),
           "Infeasible", col = "red", cex = 0.7, pos = 4, font = 3)
      text(xlim[1] + 0.4*diff(xlim), ylim[1] + 0.1*diff(ylim),
           "Feasible", col = "steelblue", cex = 0.7, pos = 4, font = 3)
    }
    
    # Add legend
    legend("topright", 
           legend = c("Start", "Bracketing", "Bisection", "Converged"),
           pch = c(16, 16, 16, 16),
           col = c("red", "pink", "gray50", "darkblue"),
           bty = "n", cex = 0.8)
    
    # Add takeaway box in a rectangle on the right side
    # Calculate position for the box - making it larger and better positioned
    box_x_left <- xlim[1] + 0.58 * diff(xlim)
    box_x_right <- xlim[1] + 0.97 * diff(xlim)
    box_y_top <- ylim[1] + 0.72 * diff(ylim)
    box_y_bottom <- ylim[1] + 0.26 * diff(ylim)
    
    # Draw rectangle
    # rect(box_x_left, box_y_bottom, box_x_right, box_y_top, 
    #      col = "white", border = "gray30", lwd = 1.5)
    
    # Add text inside rectangle with careful positioning
    takeaway_lines <- c(
      sprintf("Solution found after %d bisections,", sum(history$phase == "Bisect")),
      # "the algorithm finds a solution.",
      sprintf("Optimal t: %.5f", history$t_vals[final_idx]),
      sprintf("Final Sample Size: %s", format(history$n_vals[final_idx], big.mark = ",")),
      sprintf("(within the %s limit)", format(target, big.mark = ","))
    )
    
    # Calculate vertical spacing to fit all lines comfortably
    box_height <- box_y_top - box_y_bottom
    top_margin <- 0.10 * box_height
    bottom_margin <- 0.10 * box_height
    available_height <- box_height - top_margin - bottom_margin
    line_spacing <- available_height / (length(takeaway_lines) * 2)
    # line_spacing <- 10
    
    # Position text with left margin - using adj instead of pos for better control
    left_margin <- 0.04 * diff(xlim)
    text_x <- box_x_left + left_margin
    
    for (i in seq_along(takeaway_lines)) {
      text_y <- box_y_top - top_margin - (i - 1) * line_spacing
      text(text_x, text_y, takeaway_lines[i], 
           adj = 0, cex = 0.58, col = "black")
    }
  }
  
  if (isTRUE(plot)) {
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    
    plot_pdf_path <- if (isTRUE(save_pdf)) {
      file.path(plot_dir,
                sprintf("%s_CV_and_alloc_nmax_%d_tol_%d_t_%0.12f.pdf",
                        plot_prefix, target_size, tolerance, best_t))
    } else NULL
    
    # Base filename for PNG files
    png_base <- sprintf("%s_nmax_%d_tol_%d_t_%0.6f",
                        plot_prefix, target_size, tolerance, best_t)
    
    # Initial targets DOM+V
    targ_for_plot <- NULL
    if (!is.null(expectedCV) && isTRUE(show_targets)) {
      tmp <- errors_initial
      names(tmp) <- sub("^CV", "V", names(tmp))
      targ_for_plot <- tmp[, c("DOM", vcols_expected), drop = FALSE]
    }
    
    # Caps DOM+V (if requested)
    caps_for_plot <- NULL
    if (!is.null(expectedCV) && isTRUE(show_caps) && !is.null(caps_mat)) {
      caps_for_plot <- data.frame(DOM = rownames(caps_mat), stringsAsFactors = FALSE)
      for (k in 1:n_var) caps_for_plot[[paste0("V", k)]] <- caps_mat[, paste0("CV", k)]
    }
    
    # --- PDF ---
    if (!is.null(plot_pdf_path)) {
      pdf(plot_pdf_path, width = 11, height = 7)
      
      # NEW: Convergence plot first
      if (isTRUE(plot_convergence) && length(convergence_history$t_vals) > 0) {
        plot_convergence_trajectory(convergence_history, target_size, tolerance)
      }
      
      if (!is.null(expectedCV)) {
        for (d in expectedCV$DOM) {
          cv_vals <- as.numeric(expectedCV[expectedCV$DOM == d, vcols_expected, drop = TRUE])
          names(cv_vals) <- vcols_expected
          
          tar_vals <- NULL
          if (!is.null(targ_for_plot)) {
            tar_vals <- as.numeric(targ_for_plot[targ_for_plot$DOM == d, vcols_expected, drop = TRUE])
          }
          
          cap_vals <- NULL
          if (!is.null(caps_for_plot)) {
            cap_vals <- as.numeric(caps_for_plot[caps_for_plot$DOM == d, vcols_expected, drop = TRUE])
          }
          
          plot_cv_domain(d, cv_vals, tar_vals, cap_vals)
        }
      } else {
        plot.new()
        text(0.5, 0.5, "CV plots skipped: beat1cv_fun not provided.", cex = 1.1)
      }
      
      plot_alloc(alloc_df)
      dev.off()
    }
    
    # --- DISPLAY ---
    # NEW: Show convergence plot first
    if (isTRUE(plot_convergence) && length(convergence_history$t_vals) > 0) {
      # Save PNG
      png_file <- file.path(plot_dir, paste0(png_base, "_convergence.png"))
      png(png_file, width = 1200, height = 700, res = 120)
      plot_convergence_trajectory(convergence_history, target_size, tolerance)
      dev.off()
      message("PNG saved: ", normalizePath(png_file))
      
      # Display
      plot_convergence_trajectory(convergence_history, target_size, tolerance)
    }
    
    if (!is.null(expectedCV)) {
      for (d in expectedCV$DOM) {
        cv_vals <- as.numeric(expectedCV[expectedCV$DOM == d, vcols_expected, drop = TRUE])
        names(cv_vals) <- vcols_expected
        
        tar_vals <- NULL
        if (!is.null(targ_for_plot)) {
          tar_vals <- as.numeric(targ_for_plot[targ_for_plot$DOM == d, vcols_expected, drop = TRUE])
        }
        
        cap_vals <- NULL
        if (!is.null(caps_for_plot)) {
          cap_vals <- as.numeric(caps_for_plot[caps_for_plot$DOM == d, vcols_expected, drop = TRUE])
        }
        
        # Save PNG
        png_file <- file.path(plot_dir, sprintf("%s_CV_domain_%s.png", png_base, d))
        png(png_file, width = 1200, height = 700, res = 120)
        plot_cv_domain(d, cv_vals, tar_vals, cap_vals)
        dev.off()
        
        # Display
        plot_cv_domain(d, cv_vals, tar_vals, cap_vals)
      }
    } else {
      message("Plot CV skipped: beat1cv_fun not given.")
    }
    
    # Save allocation plot as PNG
    png_file <- file.path(plot_dir, paste0(png_base, "_allocation.png"))
    png(png_file, width = 1200, height = 700, res = 120)
    plot_alloc(alloc_df)
    dev.off()
    message("PNG saved: ", normalizePath(png_file))
    
    # Display allocation plot
    plot_alloc(alloc_df)
    
    if (!is.null(plot_pdf_path)) message("PDF saved: ", normalizePath(plot_pdf_path))
  }
  
  list(
    feasible = TRUE,
    t = best_t,
    n_total = res$n_tot,
    early_stop = early_stop,
    errors_scaled = res$errors,
    errors_initial = errors_initial,
    errors_base = errors_base,        # initial + caps applied (min)
    cv_caps = cv_caps,
    allocation = alloc_df,
    expectedCV = expectedCV,
    expectedCV_spec = expectedCV_spec,
    plot_pdf = plot_pdf_path,
    convergence_history = convergence_history  # NEW: return tracking data
  )
}
