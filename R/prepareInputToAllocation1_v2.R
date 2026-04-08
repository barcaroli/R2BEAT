# prepareInputToAllocation1_v2
#
# Extended version of prepareInputToAllocation1 that accepts domain_var as a
# character vector of ANY length, producing one DOM column per domain level
# plus the mandatory DOM1 = 1 (national constant).
#
# With domain_var = "region"  (single element, n = 1):
#   DOM1 = 1 (national)          ← same as original, fully backward-compatible
#   DOM2 = region codes
#
# With domain_var = c("region", "province")  (n = 2):
#   DOM1 = 1 (national)
#   DOM2 = region codes          ← coarsest supplied level
#   DOM3 = province codes        ← finest supplied level
#
# With domain_var = c("macro", "region", "province")  (n = 3):
#   DOM1 = 1, DOM2 = macro, DOM3 = region, DOM4 = province
#
# Key design decision
# -------------------
# buildFrameDF and buildStrataDF are called with the FINEST domain variable
# (domain_var[n]). This gives strata statistics at the finest available
# level. Coarser domain values are recovered from samp_frame via a single
# merge (valid only if domain variables are nested, which is the standard
# assumption in stratified survey design).
#
# The errors CV table passed to beat.2st must have (n + 1) rows:
#   row 1 → DOM1 national constraint
#   row 2 → DOM2 constraint (coarsest supplied domain)
#   ...
#   row n+1 → DOM(n+1) constraint (finest supplied domain)

prepareInputToAllocation1 <- function(samp_frame, id_PSU, id_SSU,
                                          strata_var, target_vars,
                                          deff_var, domain_var,
                                          minimum, delta, deff_sugg)
{
  # buildFrameDF and buildStrataDF are internal (non-exported) R2BEAT functions.
  # When prepareInputToAllocation1_v2 is sourced as a standalone file (outside
  # the R2BEAT package namespace), they must be fetched explicitly.
  buildFrameDF  <- get("buildFrameDF",  envir = asNamespace("R2BEAT"))
  buildStrataDF <- get("buildStrataDF", envir = asNamespace("R2BEAT"))

  f <- 0.05
  if (is.null(samp_frame$one))
    samp_frame$one <- 1

  # ── Multi-level domain setup ─────────────────────────────────────────────────
  domain_var <- as.character(domain_var)
  n_dom_vars <- length(domain_var)
  finest_dom <- domain_var[n_dom_vars]   # finest (last) level drives buildStrataDF

  cat("\nCalculating strata...")
  frame  <- buildFrameDF(df = samp_frame, id = id_SSU,
                         X  = strata_var, Y = target_vars,
                         domainvalue = finest_dom)
  nvarY  <- length(grep("Y", colnames(frame)))
  strata <- buildStrataDF(frame, progress = FALSE)

  # ── Build DOM columns ────────────────────────────────────────────────────────
  # After buildStrataDF, strata$DOM1 holds finest domain values.
  # Remap to:
  #   DOM1       = 1                  (national constant, always)
  #   DOM2       = domain_var[1]      (coarsest supplied, e.g. region)
  #   DOM3       = domain_var[2]      (next finer, e.g. province)
  #   ...
  #   DOM(n+1)   = domain_var[n]      (= finest_dom)

  if (n_dom_vars == 1) {
    # Original two-level behaviour — fully backward compatible
    strata$DOM2 <- strata$DOM1
    strata$DOM1 <- 1L

  } else {
    # Retrieve mapping from finest domain to all coarser domains in one merge.
    # Requires that domain variables are nested (standard survey assumption).
    coarser_vars   <- domain_var[-n_dom_vars]          # all levels except finest
    mapping_cols   <- c(finest_dom, coarser_vars)
    domain_mapping <- unique(samp_frame[, mapping_cols, drop = FALSE])

    # Ensure types match for the merge key
    strata$DOM1                    <- as.numeric(as.character(strata$DOM1))
    domain_mapping[[finest_dom]]   <- as.numeric(domain_mapping[[finest_dom]])

    # Merge: strata$DOM1 (finest domain codes) ↔ domain_mapping[finest_dom]
    strata <- merge(strata, domain_mapping,
                    by.x  = "DOM1",
                    by.y  = finest_dom,
                    all.x = TRUE)

    # DOM(n+1) = finest domain (was strata$DOM1 before merge)
    strata[[paste0("DOM", n_dom_vars + 1)]] <- strata$DOM1

    # DOM(k+1) = coarser levels in order  (k = 1 is coarsest, k = n-1 is next-finest)
    for (j_dom in seq_len(n_dom_vars - 1)) {
      strata[[paste0("DOM", j_dom + 1)]] <- strata[[ coarser_vars[j_dom] ]]
      strata[[ coarser_vars[j_dom] ]]    <- NULL     # drop the temporary column
    }

    strata$DOM1 <- 1L   # national constant
  }

  # ── Rest of function identical to original ───────────────────────────────────
  strata$STRATUM <- as.factor(strata$STRATO)
  strata$STRATO  <- strata$X1 <- NULL
  strata <- strata[order(as.numeric(as.character(strata$STRATUM))), ]

  deff <- strata$STRATUM
  deff <- as.data.frame(deff)
  colnames(deff) <- "STRATUM"
  for (i in c(1:nvarY)) {
    st <- paste0("deff$DEFF", i, " <- ", deff_sugg)
    eval(parse(text = st))
  }

  st <- paste0("b_nar <- aggregate(one ~ ", strata_var, " + ",
               id_PSU, ", data=samp_frame, FUN=sum)")
  eval(parse(text = st))
  st <- paste0("b_nar <- aggregate(one ~ ", strata_var, ", b_nar, FUN=mean)")
  eval(parse(text = st))
  if (ncol(b_nar) > 1)
    colnames(b_nar)[1] <- "STRATUM"
  if (ncol(b_nar) == 1)
    b_nar$STRATUM <- 1
  b_nar$one <- b_nar$one * f
  deff <- merge(deff, b_nar, by = "STRATUM")
  colnames(deff)[ncol(deff)] <- "b_nar"
  deff <- deff[order(as.numeric(as.character(deff$STRATUM))), ]

  effst <- strata$STRATUM
  effst <- as.data.frame(effst)
  colnames(effst)[1] <- "STRATUM"
  for (i in c(1:nvarY)) {
    st <- paste0("effst$EFFST", i, " <- 1")
    eval(parse(text = st))
  }
  effst <- effst[order(as.numeric(as.character(effst$STRATUM))), ]

  cat("\nCalculating rho in strata...")
  rho <- NULL
  rho$STRATUM <- strata$STRATUM
  rho <- as.data.frame(rho)
  for (i in c(1:nvarY)) {
    eval(parse(text = paste0("rho$RHO_AR",  i, " <- 1")))
    eval(parse(text = paste0("rho$RHO_NAR", i, " <- NA")))
  }
  L <- NULL
  k <- 0
  for (s in c(rho$STRATUM)) {
    cat("\nStratum ", s)
    k <- k + 1
    for (i in c(1:nvarY)) {
      eval(parse(text = paste0("mu <- mean(samp_frame$",
                               target_vars[i], "[samp_frame$", strata_var,
                               " == ", s, "])")))
      eval(parse(text = paste0("D2y <- sum((samp_frame$",
                               target_vars[i], "[samp_frame$", strata_var,
                               " == ", s, "]-mu)^2)")))
      eval(parse(text = paste0("L <- unique(samp_frame$",
                               id_PSU, "[samp_frame$", strata_var, "==s])")))
      D2w <- 0
      for (l in c(L)) {
        eval(parse(text = paste0("mu_L <- mean(samp_frame$",
                                 target_vars[i], "[samp_frame$", id_PSU,
                                 " == ", l, "])")))
        eval(parse(text = paste0("D2w <- D2w + sum((samp_frame$",
                                 target_vars[i], "[samp_frame$", id_PSU,
                                 " == ", l, "] - mu_L)^2)")))
      }
      eval(parse(text = paste0("rho$RHO_NAR", i, "[", k,
                               "] <- 1 - (length(L) / (length(L) - 1) * (D2w / D2y))")))
      eval(parse(text = paste0("rho$RHO_NAR", i, "[", k,
                               "] <- 1 - (D2w / D2y)")))
    }
  }

  strat_mun <- samp_frame[, c(strata_var, id_PSU)]
  strat_mun <- strat_mun[!duplicated(strat_mun), ]
  st <- paste0("mun <- aggregate(samp_frame$one,by=list(samp_frame$",
               id_PSU, "),sum)")
  eval(parse(text = st))
  colnames(mun) <- c(id_PSU, "N")
  mun      <- merge(mun, strat_mun)
  psu_file <- mun[, c(id_PSU, strata_var, "N")]
  colnames(psu_file) <- c("PSU_ID", "STRATUM", "PSU_MOS")

  des_file <- aggregate(psu_file$PSU_MOS, by = list(psu_file$STRATUM), sum)
  colnames(des_file) <- c("STRATUM", "STRAT_MOS")
  des_file$DELTA   <- delta
  des_file$MINIMUM <- minimum

  for (i in c(1:nvarY)) {
    st <- paste0("rho$RHO_NAR", i,
                 " <- ifelse(is.nan(rho$RHO_NAR", i, "),0,rho$RHO_NAR", i, ")")
    eval(parse(text = st))
  }

  out <- list(strata   = strata,
              deff     = deff,
              effst    = effst,
              rho      = rho,
              psu_file = psu_file,
              des_file = des_file)
  return(out)
}
