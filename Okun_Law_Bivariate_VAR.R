# =============================================================================
# Output Growth and Unemployment Dynamics in the United States
# A Bivariate VAR Test of Okun's Law, 1984Q1-2019Q4
#
# Author : Sameer Chawla (MSc Economics, GIPE Pune)
# Data   : FRED - GDPC1 (real GDP, quarterly), UNRATE (unemployment rate, monthly)
#
# Pipeline
#   0. Setup and data            5. Granger and instantaneous causality
#   1. Descriptive statistics    6. Orthogonalised IRFs
#   2. Unit-root tests           7. FEVD
#   3. Lag-order audit           8. Robustness
#   4. Baseline VAR(3) + checks  9. Static Okun benchmark (OLS + Newey-West)
#
# Outputs: all figures (.png) and tables (.csv) are written to ./output/
# NOTE: model calls use a LITERAL p (p = 3L). This is deliberate: vars::irf()
# with boot = TRUE re-evaluates the stored VAR call, and a p passed through a
# variable or loop index can fail to resolve. To change the baseline lag,
# edit every `p = 3L` in sections 4-8.
# =============================================================================

# ---- 0. SETUP ----------------------------------------------------------------
MAX_P       <- 6L      # maximum lag considered in lag-order selection
H           <- 12L     # IRF / FEVD horizon in quarters
GRANGER_BOOT <- 4999L  # bootstrap replications, Granger tests (use 199 to test)
IRF_RUNS    <- 1999L   # bootstrap replications, IRF bands (use 199 to test)
SEED        <- 20261020L
OUT         <- "output"

needed  <- c("quantmod", "zoo", "urca", "vars", "sandwich", "lmtest")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(quantmod); library(zoo); library(urca)
  library(vars); library(sandwich); library(lmtest)
})

dir.create(OUT, showWarnings = FALSE)
options(width = 115)
set.seed(SEED)

heading <- function(x) cat("\n\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
pval    <- function(x) as.numeric(x$p.value)
fmt     <- function(x, d = 4) round(x, d)
save_csv <- function(x, file) write.csv(x, file.path(OUT, file), row.names = FALSE)
save_png <- function(file, expr, width = 1800, height = 1100, res = 200) {
  png(file.path(OUT, file), width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  invisible(NULL)
}

heading("0. SOFTWARE AND DATA SOURCE")
cat("Data retrieved on:", format(Sys.Date()), "\n")
cat("R version:", R.version.string, "\n")
print(sapply(needed, function(z) as.character(packageVersion(z))))
cat("FRED series: GDPC1 (quarterly real GDP); UNRATE (monthly unemployment rate).\n")
cat("Sample: 1984Q1 to 2019Q4; monthly UNRATE aggregated to quarterly mean.\n")
cat("GDP growth = 100 x log-difference (percent per quarter); dU = change in rate (pp).\n")

fred_env <- new.env(parent = globalenv())
quantmod::getSymbols(c("GDPC1", "UNRATE"), src = "FRED", env = fred_env, auto.assign = TRUE)

gdpq <- zoo::zoo(as.numeric(fred_env$GDPC1),
                 order.by = zoo::as.yearqtr(zoo::index(fred_env$GDPC1)))
uneq <- stats::aggregate(
  zoo::zoo(as.numeric(fred_env$UNRATE), order.by = as.Date(zoo::index(fred_env$UNRATE))),
  by = zoo::as.yearqtr, FUN = mean
)
lev <- na.omit(merge(GDP = gdpq, UNRATE = uneq, all = FALSE))
lev <- window(lev, start = zoo::as.yearqtr("1983 Q4"), end = zoo::as.yearqtr("2019 Q4"))
stopifnot(nrow(lev) == 145L, !anyNA(lev))

g   <- 100 * diff(log(lev[, "GDP"]))
du  <- diff(lev[, "UNRATE"])
dat <- na.omit(merge(GDP_Growth = g, Delta_Unemployment = du))
stopifnot(nrow(dat) == 144L, as.character(zoo::index(dat)[1]) == "1984 Q1")

Y <- ts(zoo::coredata(dat), start = c(1984, 1), frequency = 4)
colnames(Y) <- c("GDP_Growth", "Delta_Unemployment")

# ---- 1. DESCRIPTIVE STATISTICS ------------------------------------------------
heading("1. DESCRIPTIVE STATISTICS")
desc <- data.frame(variable = colnames(Y), n = nrow(Y), mean = apply(Y, 2, mean),
                   sd = apply(Y, 2, sd), min = apply(Y, 2, min), max = apply(Y, 2, max))
print(desc, row.names = FALSE)
save_csv(desc, "table_descriptives.csv")
cat("\nContemporaneous Pearson correlation:\n"); print(cor(Y))
save_png("01_series.png",
         plot(Y, main = "US GDP growth and unemployment-rate change (quarterly)", xlab = "Year"))

# ---- 2. UNIT-ROOT TESTS: ADF, PP, KPSS ---------------------------------------
heading("2. UNIT-ROOT TESTS: ADF, PP (null = unit root) AND KPSS (null = stationarity)")
results <- list()
for (v in colnames(Y)) {
  x <- as.numeric(Y[, v])
  cat("\n=======================================\nVARIABLE:", v, "\n=======================================\n")
  adf      <- urca::ur.df(x, type = "drift", lags = 8, selectlags = "AIC")
  pp_short <- urca::ur.pp(x, type = "Z-tau", model = "constant", lags = "short")
  pp_long  <- urca::ur.pp(x, type = "Z-tau", model = "constant", lags = "long")
  kp_short <- urca::ur.kpss(x, type = "mu", lags = "short")
  kp_long  <- urca::ur.kpss(x, type = "mu", lags = "long")
  cat("\n----- ADF -----\n");        print(summary(adf))
  cat("\n----- PP short -----\n");   print(summary(pp_short))
  cat("\n----- PP long -----\n");    print(summary(pp_long))
  cat("\n----- KPSS short -----\n"); print(summary(kp_short))
  cat("\n----- KPSS long -----\n");  print(summary(kp_long))
  results[[v]] <- data.frame(
    Variable = v,
    ADF = as.numeric(adf@teststat[1]),           ADF_crit5 = as.numeric(adf@cval[1, "5pct"]),
    PP_short = as.numeric(pp_short@teststat),    PP_long = as.numeric(pp_long@teststat),
    PP_crit5 = as.numeric(pp_short@cval[1, "5pct"]),
    KPSS_short = as.numeric(kp_short@teststat),  KPSS_long = as.numeric(kp_long@teststat),
    KPSS_crit5 = as.numeric(kp_long@cval[1, "5pct"])
  )
}
stationarity_table <- do.call(rbind, results)
cat("\nCOMPARATIVE STATIONARITY TABLE:\n")
print(data.frame(Variable = stationarity_table$Variable, round(stationarity_table[, -1], 4)),
      row.names = FALSE)
save_csv(stationarity_table, "table_unit_roots.csv")
cat("\nRule: ADF/PP reject a unit root if statistic < 5% critical value.\n",
    "      KPSS rejects stationarity if statistic > 5% critical value.\n", sep = "")

# ---- 3. LAG-ORDER AUDIT -------------------------------------------------------
heading("3. LAG-ORDER SELECTION AND COMMON-SAMPLE COMPARISON")
sel <- vars::VARselect(Y, lag.max = MAX_P, type = "const")
cat("Information-criterion recommendations:\n"); print(sel$selection)
cat("\nFull criteria (VARselect holds the effective sample fixed across p):\n")
print(round(sel$criteria, 4))
save_csv(data.frame(criterion = rownames(sel$criteria), round(sel$criteria, 4)),
         "table_lag_criteria.csv")

# VAR() drops p observations, so fits are aligned to the SAME 138 effective
# quarters (original rows 7:144) to compare lag orders fairly.
y1 <- Y[6:144, , drop = FALSE]; y2 <- Y[5:144, , drop = FALSE]
y3 <- Y[4:144, , drop = FALSE]; y4 <- Y[3:144, , drop = FALSE]
m1 <- vars::VAR(y1, p = 1L, type = "const"); m2 <- vars::VAR(y2, p = 2L, type = "const")
m3 <- vars::VAR(y3, p = 3L, type = "const"); m4 <- vars::VAR(y4, p = 4L, type = "const")
cs_models <- list(VAR1 = m1, VAR2 = m2, VAR3 = m3, VAR4 = m4)
stopifnot(all(vapply(cs_models, function(m) nrow(residuals(m)), numeric(1)) == 138L))

# HC1 heteroskedasticity-robust covariance built from the VAR's own design matrix
hc1_var_cov <- function(fit, response) {
  eq <- fit$varresult[[response]]
  cn <- names(coef(eq))
  X  <- as.matrix(fit$datamat[, cn, drop = FALSE]); storage.mode(X) <- "double"
  u  <- as.numeric(residuals(eq))
  n  <- nrow(X); k <- ncol(X)
  stopifnot(n == length(u), all(is.finite(X)), all(is.finite(u)), qr(X)$rank == k)
  XXi <- solve(crossprod(X))
  S   <- crossprod(X, sweep(X, 1, u^2, "*"))
  V   <- (n / (n - k)) * XXi %*% S %*% XXi
  dimnames(V) <- list(cn, cn)
  V
}
# Joint Wald test: all p lags of `cause` are zero in the `response` equation
robust_wald <- function(fit, response, cause) {
  b  <- coef(fit$varresult[[response]])
  V  <- hc1_var_cov(fit, response)
  ii <- match(paste0(cause, ".l", seq_len(fit$p)), names(b))
  if (anyNA(ii)) stop("Cross-lag coefficients not found.")
  beta <- matrix(b[ii], ncol = 1)
  W    <- as.numeric(t(beta) %*% solve(V[ii, ii, drop = FALSE], beta))
  data.frame(Chi_Square = W, DF = length(ii), P_Value = pchisq(W, df = length(ii), lower.tail = FALSE))
}

lag_table <- do.call(rbind, lapply(names(cs_models), function(nm) {
  m <- cs_models[[nm]]
  a <- robust_wald(m, "Delta_Unemployment", "GDP_Growth")
  b <- robust_wald(m, "GDP_Growth", "Delta_Unemployment")
  data.frame(
    model = nm, n = nrow(residuals(m)),
    BG_p   = pval(vars::serial.test(m, lags.bg = 8, type = "BG")$serial),
    JB_p   = pval(vars::normality.test(m, multivariate.only = TRUE)$jb.mul$JB),
    ARCH_p = pval(vars::arch.test(m, lags.multi = 2, multivariate.only = TRUE)$arch.mul),
    max_root = max(Mod(vars::roots(m))),
    classical_GDP_to_dU = pval(vars::causality(m, cause = "GDP_Growth")$Granger),
    HC1_GDP_to_dU       = a$P_Value,
    classical_dU_to_GDP = pval(vars::causality(m, cause = "Delta_Unemployment")$Granger),
    HC1_dU_to_GDP       = b$P_Value
  )
}))
cat("\nCommon-sample comparison (n = 138 for every model):\n")
print(transform(lag_table, BG_p = fmt(BG_p), JB_p = fmt(JB_p), ARCH_p = fmt(ARCH_p),
                max_root = fmt(max_root), classical_GDP_to_dU = signif(classical_GDP_to_dU, 3),
                HC1_GDP_to_dU = signif(HC1_GDP_to_dU, 3),
                classical_dU_to_GDP = signif(classical_dU_to_GDP, 3),
                HC1_dU_to_GDP = signif(HC1_dU_to_GDP, 3)), row.names = FALSE)
save_csv(lag_table, "table_lag_common_sample.csv")
cat("\nBaseline choice: VAR(3) (AIC and FPE optimum; parsimonious vs VAR(4)).\n",
    "VAR(4) is retained as a robustness check. Failure to reject a diagnostic\n",
    "null does not prove correct specification.\n", sep = "")

# ---- 4. BASELINE VAR(3): ESTIMATION AND DIAGNOSTICS ---------------------------
heading("4. BASELINE VAR(3): ESTIMATES AND DIAGNOSTICS (FULL SAMPLE, 141 OBS)")
model <- vars::VAR(Y, p = 3L, type = "const")
print(summary(model))

cat("\nBreusch-Godfrey LM (8 lags), H0: no residual serial correlation:\n")
print(vars::serial.test(model, lags.bg = 8, type = "BG")$serial)
cat("\nPortmanteau (16 lags), H0: no residual serial correlation:\n")
print(vars::serial.test(model, lags.pt = 16, type = "PT.asymptotic")$serial)
cat("\nMultivariate Jarque-Bera, H0: normality:\n")
print(vars::normality.test(model, multivariate.only = TRUE)$jb.mul)
cat("\nMultivariate ARCH-LM (2 lags), H0: no ARCH:\n")
print(vars::arch.test(model, lags.multi = 2, multivariate.only = TRUE)$arch.mul)
cat("\nCompanion eigenvalue moduli (stable if all < 1):\n")
print(Mod(vars::roots(model)))
save_png("02_cusum.png", plot(vars::stability(model, type = "OLS-CUSUM")))

diag_baseline <- data.frame(
  test = c("Breusch-Godfrey (8)", "Jarque-Bera", "ARCH-LM (2)"),
  p_value = c(pval(vars::serial.test(model, lags.bg = 8, type = "BG")$serial),
              pval(vars::normality.test(model, multivariate.only = TRUE)$jb.mul$JB),
              pval(vars::arch.test(model, lags.multi = 2, multivariate.only = TRUE)$arch.mul)))
save_csv(diag_baseline, "table_diagnostics_baseline.csv")
cat("\nResidual correlation:\n"); print(cor(residuals(model)))
# Save the VAR residual scatter shown in the README / appendix.
save_png("10_var_residuals.png", {
  U <- residuals(model)
  plot(U[, "GDP_Growth"], U[, "Delta_Unemployment"],
       xlab = "GDP growth equation residual",
       ylab = "Unemployment-change equation residual (pp)",
       main = "Reduced-form VAR(3) residuals",
       pch = 1, cex = 0.85)
  abline(h = 0, v = 0, col = "gray80", lty = 3)
})

# ---- 5. GRANGER AND INSTANTANEOUS CAUSALITY -----------------------------------
heading("5. GRANGER PREDICTABILITY AND INSTANTANEOUS CAUSALITY (VAR(3))")
cat("Null: all cross-lag coefficients are zero. Three p-values per direction:\n",
    "classical F, bootstrap F (robust to non-normality/heteroskedasticity),\n",
    "and HC1-robust Wald chi-square.\n", sep = "")
set.seed(SEED)
c_g <- vars::causality(model, cause = "GDP_Growth", boot = TRUE, boot.runs = GRANGER_BOOT)
set.seed(SEED + 1L)
c_u <- vars::causality(model, cause = "Delta_Unemployment", boot = TRUE, boot.runs = GRANGER_BOOT)
c_g_cl <- vars::causality(model, cause = "GDP_Growth")
c_u_cl <- vars::causality(model, cause = "Delta_Unemployment")
w_g <- robust_wald(model, "Delta_Unemployment", "GDP_Growth")
w_u <- robust_wald(model, "GDP_Growth", "Delta_Unemployment")

cat("\nGDP -> dU (bootstrap):\n");  print(c_g$Granger)
cat("\ndU -> GDP (bootstrap):\n");  print(c_u$Granger)
cat("\nInstantaneous causality (residual covariance test):\n"); print(c_g$Instant)

granger_table <- data.frame(
  direction = c("GDP -> dU", "dU -> GDP"),
  F_stat = fmt(c(as.numeric(c_g$Granger$statistic), as.numeric(c_u$Granger$statistic))),
  classical_p = signif(c(pval(c_g_cl$Granger), pval(c_u_cl$Granger)), 3),
  bootstrap_p = fmt(c(pval(c_g$Granger), pval(c_u$Granger))),
  HC1_Wald_chisq = fmt(c(w_g$Chi_Square, w_u$Chi_Square)),
  HC1_Wald_p = signif(c(w_g$P_Value, w_u$P_Value), 3)
)
cat("\nGranger summary:\n"); print(granger_table, row.names = FALSE)
save_csv(granger_table, "table_granger.csv")
cat("\nGranger causality is predictive precedence, not structural causation.\n")

# ---- 6. ORTHOGONALISED IRFs ---------------------------------------------------
heading("6. ORTHOGONALISED IRFs, 95% BOOTSTRAP BANDS (ORDER: GDP FIRST)")
cat("Impact response of dU to a GDP shock may be non-zero; the reverse impact\n",
    "response is zero by construction under this ordering.\n", sep = "")
irf_g <- vars::irf(model, impulse = "GDP_Growth", response = "Delta_Unemployment",
                   n.ahead = H, ortho = TRUE, boot = TRUE, runs = IRF_RUNS, ci = .95,
                   seed = SEED + 2L)
save_png("03_irf_gdp_to_unemployment.png",
         plot(irf_g, main = "Response of change in unemployment to GDP growth shock"))
g_irf <- data.frame(horizon = 0:H, estimate = irf_g$irf$GDP_Growth[, 1],
                    lower = irf_g$Lower$GDP_Growth[, 1], upper = irf_g$Upper$GDP_Growth[, 1])
g_irf$CI_excludes_zero <- with(g_irf, lower > 0 | upper < 0)
cat("\nGDP shock -> dU:\n"); print(transform(g_irf, estimate = fmt(estimate), lower = fmt(lower),
                                              upper = fmt(upper)), row.names = FALSE)
save_csv(g_irf, "table_irf_gdp_to_unemployment.csv")

irf_u <- vars::irf(model, impulse = "Delta_Unemployment", response = "GDP_Growth",
                   n.ahead = H, ortho = TRUE, boot = TRUE, runs = IRF_RUNS, ci = .95,
                   seed = SEED + 3L)
save_png("04_irf_unemployment_to_gdp.png",
         plot(irf_u, main = "Response of GDP growth to unemployment-change shock"))
u_irf <- data.frame(horizon = 0:H, estimate = irf_u$irf$Delta_Unemployment[, 1],
                    lower = irf_u$Lower$Delta_Unemployment[, 1],
                    upper = irf_u$Upper$Delta_Unemployment[, 1])
cat("\ndU shock -> GDP:\n"); print(transform(u_irf, estimate = fmt(estimate), lower = fmt(lower),
                                              upper = fmt(upper)), row.names = FALSE)
save_csv(u_irf, "table_irf_unemployment_to_gdp.csv")

# ---- 7. FEVD ------------------------------------------------------------------
heading("7. FORECAST ERROR VARIANCE DECOMPOSITION (VAR(3), GDP FIRST)")
cat("Shares are fractions of forecast-error variance, not of realised outcomes.\n")
fe <- vars::fevd(model, n.ahead = H)
for (v in names(fe)) {
  cat("\nFEVD for ", v, " (%):\n", sep = ""); print(round(100 * fe[[v]], 2))
}
# Publication layout: two stacked FEVD panels and a separate legend strip.
# The bar heights are percentages of FORECAST-ERROR VARIANCE.
plot_fevd_clean <- function(fe) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  layout(matrix(c(1, 2, 3), ncol = 1), heights = c(1, 1, 0.22))
  for (v in c("GDP_Growth", "Delta_Unemployment")) {
    share_matrix <- t(100 * fe[[v]][, c("GDP_Growth", "Delta_Unemployment"), drop = FALSE])
    par(mar = c(4.1, 4.8, 3.2, 1.0))
    barplot(share_matrix,
            col = c("gray25", "gray82"), border = "white",
            space = 0.25, ylim = c(0, 100),
            names.arg = seq_len(ncol(share_matrix)),
            xlab = "Forecast horizon (quarters)",
            ylab = "FE variance share (%)",
            main = paste("FEVD for", v), cex.main = 1.05)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = c("GDP growth innovation", "Unemployment-change innovation"),
         fill = c("gray25", "gray82"), horiz = TRUE, bty = "n", cex = 0.9)
  invisible(NULL)
}
save_png("05_fevd.png", plot_fevd_clean(fe), width = 1900, height = 1500, res = 190)
save_csv(data.frame(horizon = 1:H, round(100 * fe$Delta_Unemployment, 2)), "table_fevd_unemployment.csv")
save_csv(data.frame(horizon = 1:H, round(100 * fe$GDP_Growth, 2)), "table_fevd_gdp.csv")

# ---- 8. ROBUSTNESS ------------------------------------------------------------
heading("8A. REVERSE CHOLESKY ORDERING (DU FIRST)")
Y_rev <- Y[, c("Delta_Unemployment", "GDP_Growth")]
model_rev <- vars::VAR(Y_rev, p = 3L, type = "const")
irf_rev <- vars::irf(model_rev, impulse = "GDP_Growth", response = "Delta_Unemployment",
                     n.ahead = H, ortho = TRUE, boot = TRUE, runs = IRF_RUNS, ci = .95,
                     seed = SEED + 3L)
save_png("06_irf_reverse_ordering.png",
         plot(irf_rev, main = "GDP shock -> unemployment: reverse Cholesky ordering"))
fe_rev <- vars::fevd(model_rev, n.ahead = H)
order_table <- data.frame(
  horizon = 0:H,
  IRF_GDPfirst = fmt(irf_g$irf$GDP_Growth[, 1]),
  IRF_dUfirst  = fmt(irf_rev$irf$GDP_Growth[, 1]),
  FEVD_GDPfirst = c(NA, round(100 * fe$Delta_Unemployment[, "GDP_Growth"], 2)),
  FEVD_dUfirst  = c(NA, round(100 * fe_rev$Delta_Unemployment[, "GDP_Growth"], 2))
)
print(order_table, row.names = FALSE)
save_csv(order_table, "table_ordering_robustness.csv")
cat("\nThe two orderings define DIFFERENT orthogonalised shocks. The Granger result\n",
    "does not depend on ordering; the size of the IRF and the FEVD share do.\n", sep = "")

heading("8B. LAG SENSITIVITY OF THE GDP -> dU RESPONSE (COMMON SAMPLE)")
ig1 <- vars::irf(m1, impulse = "GDP_Growth", response = "Delta_Unemployment", n.ahead = 8L,
                 ortho = TRUE, boot = TRUE, runs = IRF_RUNS, ci = .95, seed = SEED + 5L)
ig3 <- vars::irf(m3, impulse = "GDP_Growth", response = "Delta_Unemployment", n.ahead = 8L,
                 ortho = TRUE, boot = TRUE, runs = IRF_RUNS, ci = .95, seed = SEED + 6L)
ig4 <- vars::irf(m4, impulse = "GDP_Growth", response = "Delta_Unemployment", n.ahead = 8L,
                 ortho = TRUE, boot = TRUE, runs = IRF_RUNS, ci = .95, seed = SEED + 7L)
lag_irf <- data.frame(horizon = 0:8,
  VAR1 = fmt(ig1$irf$GDP_Growth[, 1]), VAR3 = fmt(ig3$irf$GDP_Growth[, 1]),
  VAR4 = fmt(ig4$irf$GDP_Growth[, 1]),
  VAR1_excl0 = with(list(l = ig1$Lower$GDP_Growth[, 1], u = ig1$Upper$GDP_Growth[, 1]), l > 0 | u < 0),
  VAR3_excl0 = with(list(l = ig3$Lower$GDP_Growth[, 1], u = ig3$Upper$GDP_Growth[, 1]), l > 0 | u < 0),
  VAR4_excl0 = with(list(l = ig4$Lower$GDP_Growth[, 1], u = ig4$Upper$GDP_Growth[, 1]), l > 0 | u < 0))
print(lag_irf, row.names = FALSE)
save_csv(lag_irf, "table_lag_sensitivity_irf.csv")
save_png("09_lag_sensitivity_irf.png", {
  rng <- range(ig1$irf$GDP_Growth[, 1], ig3$irf$GDP_Growth[, 1], ig4$irf$GDP_Growth[, 1])
  plot(0:8, ig1$irf$GDP_Growth[, 1], type = "l", lwd = 2, lty = 2, ylim = rng,
       xlab = "Quarter horizon (0 = impact)", ylab = "Response of change in unemployment (pp)",
       main = "GDP shock -> unemployment: common-sample lag sensitivity")
  lines(0:8, ig3$irf$GDP_Growth[, 1], lwd = 2, lty = 3)
  lines(0:8, ig4$irf$GDP_Growth[, 1], lwd = 2, lty = 1)
  abline(h = 0, col = "gray55")
  legend("bottomright", legend = c("VAR(1)", "VAR(3)", "VAR(4)"), lty = c(2, 3, 1), lwd = 2, bty = "n")
})

heading("8C. PRE-GFC SUBSAMPLE (THROUGH 2007Q4), VAR(3)")
Y_pre <- window(Y, end = c(2007, 4))
model_pre <- vars::VAR(Y_pre, p = 3L, type = "const")
cat("Observations:", nrow(Y_pre), "| effective:", nrow(residuals(model_pre)), "\n")
pre_diag <- data.frame(
  BG   = pval(vars::serial.test(model_pre, lags.bg = 8, type = "BG")$serial),
  JB   = pval(vars::normality.test(model_pre, multivariate.only = TRUE)$jb.mul$JB),
  ARCH = pval(vars::arch.test(model_pre, lags.multi = 2, multivariate.only = TRUE)$arch.mul))
print(pre_diag)
cat("\nPre-GFC Granger tests (classical F):\n")
print(vars::causality(model_pre, cause = "GDP_Growth")$Granger)
print(vars::causality(model_pre, cause = "Delta_Unemployment")$Granger)
irf_pre <- vars::irf(model_pre, impulse = "GDP_Growth", response = "Delta_Unemployment",
                     n.ahead = H, ortho = TRUE, boot = TRUE, runs = IRF_RUNS, ci = .95,
                     seed = SEED + 8L)
save_png("07_irf_pre_gfc.png",
         plot(irf_pre, main = "GDP shock -> unemployment: pre-GFC subsample"))
pre_irf <- data.frame(horizon = 0:H, estimate = irf_pre$irf$GDP_Growth[, 1],
                      lower = irf_pre$Lower$GDP_Growth[, 1], upper = irf_pre$Upper$GDP_Growth[, 1])
print(transform(pre_irf, estimate = fmt(estimate), lower = fmt(lower), upper = fmt(upper)), row.names = FALSE)
save_csv(pre_irf, "table_irf_pre_gfc.csv")

heading("8D. GFC IMPULSE DUMMIES (2008Q4, 2009Q1), VAR(3)")
qtime <- time(Y)
D <- cbind(D2008Q4 = as.integer(abs(qtime - 2008.75) < 1e-8),
           D2009Q1 = as.integer(abs(qtime - 2009) < 1e-8))
stopifnot(sum(D) == 2)
model_dum <- vars::VAR(Y, p = 3L, type = "const", exogen = D)
dum_diag <- data.frame(
  BG   = pval(vars::serial.test(model_dum, lags.bg = 8, type = "BG")$serial),
  JB   = pval(vars::normality.test(model_dum, multivariate.only = TRUE)$jb.mul$JB),
  ARCH = pval(vars::arch.test(model_dum, lags.multi = 2, multivariate.only = TRUE)$arch.mul))
print(dum_diag)
cat("Impulse dummies are an intervention sensitivity test, not proof of correct specification.\n")

robust_diag <- rbind(
  data.frame(spec = "Baseline VAR(3), full sample", BG = diag_baseline$p_value[1],
             JB = diag_baseline$p_value[2], ARCH = diag_baseline$p_value[3]),
  data.frame(spec = "Pre-GFC (to 2007Q4)", pre_diag),
  data.frame(spec = "GFC impulse dummies", dum_diag))
cat("\nDiagnostic p-values across specifications:\n")
print(transform(robust_diag, BG = fmt(BG), JB = fmt(JB), ARCH = fmt(ARCH)), row.names = FALSE)
save_csv(robust_diag, "table_robustness_diagnostics.csv")

# ---- 9. STATIC OKUN BENCHMARK -------------------------------------------------
heading("9. STATIC DIFFERENCE-FORM OKUN LAW: OLS WITH NEWEY-WEST STANDARD ERRORS")
df <- data.frame(g = as.numeric(Y[, "GDP_Growth"]), du = as.numeric(Y[, "Delta_Unemployment"]))
okun <- lm(du ~ g, data = df)
cat("Equation: dU = alpha + beta * g + error\n"); print(summary(okun))
Vnw <- sandwich::NeweyWest(okun, lag = 4L, prewhite = FALSE, adjust = TRUE)
cat("\nNewey-West HAC (4 lags, no prewhitening):\n"); print(lmtest::coeftest(okun, vcov. = Vnw))
ci_nw <- cbind(estimate = coef(okun),
               lower = coef(okun) - qnorm(.975) * sqrt(diag(Vnw)),
               upper = coef(okun) + qnorm(.975) * sqrt(diag(Vnw)))
cat("\nApproximate 95% HAC CI:\n"); print(round(ci_nw, 4))
g_star <- -coef(okun)[1] / coef(okun)[2]
cat("\nGrowth rate consistent with a stable unemployment rate (-alpha/beta):",
    round(g_star, 3), "% per quarter (about", round(4 * g_star, 2), "% annualised).\n")
cat("beta is a contemporaneous association, not structural causation.\n")
save_csv(data.frame(term = rownames(ci_nw), round(ci_nw, 4)), "table_okun_static_hac.csv")
save_png("08_okun_scatter.png", {
  plot(df$g, df$du, pch = 19, cex = .7, col = "gray45",
       xlab = "Quarterly GDP log growth (%)", ylab = "Quarterly change in unemployment rate (pp)",
       main = "Difference-form Okun relationship: US, 1984Q1-2019Q4")
  abline(okun, col = "firebrick", lwd = 2); abline(h = 0, v = 0, col = "gray80", lty = 3)
})

heading("DONE")
cat("Figures (.png) and tables (.csv) written to:", normalizePath(OUT), "\n")
cat("FRED data are revised over time (especially GDPC1); numbers may differ slightly on re-runs.\n")
print(sessionInfo())