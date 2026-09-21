# ============================================================================
# Output Growth and Unemployment Dynamics in the United States:
# A Bivariate VAR Test of Okun's Law, 1984Q1-2019Q4
#
# Author:  Sameer Chawla
# Program: MSc Economics, Gokhale Institute of Politics and Economics (GIPE)
# ============================================================================
#
# Pipeline:
#   1. Data construction and descriptive statistics
#   2. Unit root tests (ADF, KPSS)
#   3. Lag order selection (information criteria, then Breusch-Godfrey check)
#   4. VAR(p) estimation and residual diagnostics
#   5. Stability (characteristic roots, OLS-CUSUM)
#   6. Granger and instantaneous causality
#   7. Orthogonalised impulse responses (Cholesky, bootstrap CI)
#   8. Forecast error variance decomposition
#   9. Robustness: reverse ordering, pre-GFC subsample, GFC dummies,
#      lag sensitivity
# ============================================================================

# Install any missing packages, then load them
for (pkg in c("quantmod", "zoo", "urca", "vars")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(quantmod)
library(zoo)
library(urca)
library(vars)

set.seed(2026)
dir.create("output", showWarnings = FALSE)

# Helper: save a plot to output/
save_png <- function(name, expr, w = 1600, h = 1100, res = 170) {
  png(file.path("output", name), width = w, height = h, res = res)
  on.exit(dev.off())
  expr
}

# Settings
SAMPLE_START <- "1984 Q1"
SAMPLE_END   <- "2019 Q4"
LAG_MAX      <- 6      # maximum lag considered in Stage 3
IRF_HORIZON  <- 12     # quarters
IRF_RUNS     <- 1000   # bootstrap replications (lower for a quick test run)

# ----------------------------------------------------------------------------
# 1. Data construction
# ----------------------------------------------------------------------------
getSymbols(c("GDPC1", "UNRATE"), src = "FRED")

# Real GDP: quarterly
gdp_q <- zoo(as.numeric(GDPC1), order.by = as.yearqtr(index(GDPC1)))

# Unemployment rate: monthly -> quarterly average
unemp_q <- aggregate(zoo(as.numeric(UNRATE), order.by = index(UNRATE)),
                     as.yearqtr, mean)

macro <- merge(GDP = gdp_q, UNRATE = unemp_q, all = FALSE)

# Transformations
gdp_growth  <- 100 * diff(log(macro[, "GDP"]))   # GDP growth, % per quarter
delta_unemp <- diff(macro[, "UNRATE"])           # change in unemployment rate

var_data <- na.omit(merge(GDP_Growth = gdp_growth,
                          Delta_Unemployment = delta_unemp))
var_data <- window(var_data, start = as.yearqtr(SAMPLE_START),
                   end = as.yearqtr(SAMPLE_END))

Y <- ts(coredata(var_data), start = c(1984, 1), frequency = 4)
colnames(Y) <- c("GDP_Growth", "Delta_Unemployment")

stopifnot(nrow(Y) == 144)   # 1984Q1-2019Q4
cat("Sample:", SAMPLE_START, "-", SAMPLE_END, "| Observations:", nrow(Y), "\n")
cat("Data retrieved on", format(Sys.Date()),
    "- FRED series are revised, so results can change slightly over time.\n")

save_png("01_series.png",
         plot(Y, main = "US GDP Growth and Change in Unemployment", xlab = "Year"))

desc <- data.frame(
  Variable = colnames(Y),
  Mean = round(apply(Y, 2, mean), 3),
  SD   = round(apply(Y, 2, sd), 3),
  Min  = round(apply(Y, 2, min), 3),
  Max  = round(apply(Y, 2, max), 3)
)
print(desc, row.names = FALSE)
write.csv(desc, "output/01_descriptives.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
# 2. Unit root tests
#    ADF: H0 = unit root.   KPSS: H0 = stationary.
# ----------------------------------------------------------------------------
unit_root_row <- function(x, label) {
  adf  <- ur.df(x, type = "drift", lags = 4, selectlags = "AIC")
  kpss <- ur.kpss(x, type = "mu", lags = "long")
  data.frame(
    series = label,
    adf_stat = round(adf@teststat[1], 3),
    adf_crit_5pct = adf@cval[1, "5pct"],
    adf_rejects_unit_root = adf@teststat[1] < adf@cval[1, "5pct"],
    kpss_stat = round(kpss@teststat, 3),
    kpss_crit_5pct = kpss@cval[1, "5pct"],
    kpss_rejects_stationarity = kpss@teststat > kpss@cval[1, "5pct"]
  )
}
ur_tab <- do.call(rbind, lapply(colnames(Y), function(v) unit_root_row(Y[, v], v)))
print(ur_tab, row.names = FALSE)
write.csv(ur_tab, "output/02_unit_roots.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
# 3. Lag order selection
#    Start from the SC (BIC) choice; raise p until the Breusch-Godfrey test
#    (8 lags) no longer rejects residual serial correlation.
# ----------------------------------------------------------------------------
lag_selection <- VARselect(Y, lag.max = LAG_MAX, type = "const")
print(lag_selection$selection)
write.csv(round(lag_selection$criteria, 4), "output/03_lag_criteria.csv")

bg_pvalue <- function(fit) serial.test(fit, lags.bg = 8, type = "BG")$serial$p.value

p <- max(1, as.integer(lag_selection$selection["SC(n)"]))
var_model <- VAR(Y, p = p, type = "const")
serial_p  <- bg_pvalue(var_model)
bg_trace  <- data.frame(lag = p, bg_p_value = round(serial_p, 4))

while (serial_p < 0.05 && p < LAG_MAX) {
  p <- p + 1
  var_model <- VAR(Y, p = p, type = "const")
  serial_p  <- bg_pvalue(var_model)
  bg_trace  <- rbind(bg_trace, data.frame(lag = p, bg_p_value = round(serial_p, 4)))
}
print(bg_trace, row.names = FALSE)
write.csv(bg_trace, "output/03_lag_refinement.csv", row.names = FALSE)

if (serial_p < 0.05) {
  warning("Serial correlation remains up to LAG_MAX. Report this as a limitation.")
}
cat("Final VAR lag order p =", p, "\n")

# ----------------------------------------------------------------------------
# 4. Estimation and residual diagnostics
# ----------------------------------------------------------------------------
writeLines(capture.output(summary(var_model)), "output/04_var_estimates.txt")

sr <- serial.test(var_model, lags.bg = 8, type = "BG")
nm <- normality.test(var_model, multivariate.only = TRUE)
ar <- arch.test(var_model, lags.multi = 2, multivariate.only = TRUE)

diag_tab <- data.frame(
  test = c("Breusch-Godfrey (no serial correlation)",
           "Jarque-Bera (multivariate normality)",
           "ARCH-LM (no heteroskedasticity)"),
  p_value = round(c(sr$serial$p.value, nm$jb.mul$JB$p.value, ar$arch.mul$p.value), 4)
)
diag_tab$verdict <- ifelse(diag_tab$p_value > 0.05, "pass", "FAIL")
print(diag_tab, row.names = FALSE)
write.csv(diag_tab, "output/04_diagnostics.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
# 5. Stability
# ----------------------------------------------------------------------------
var_roots <- roots(var_model)
cat("Largest root modulus:", round(max(var_roots), 4), "\n")
if (max(var_roots) >= 1) warning("VAR is not dynamically stable.")

save_png("02_cusum.png", plot(stability(var_model, type = "OLS-CUSUM")), h = 1500)

# ----------------------------------------------------------------------------
# 6. Granger and instantaneous causality
#    Bootstrap p-values because the full-sample residuals are non-normal.
# ----------------------------------------------------------------------------
set.seed(2026)
gr_gdp   <- causality(var_model, cause = "GDP_Growth", boot = TRUE, boot.runs = 1000)
set.seed(2026)
gr_unemp <- causality(var_model, cause = "Delta_Unemployment", boot = TRUE, boot.runs = 1000)

granger_tab <- data.frame(
  hypothesis = c("GDP_Growth does not Granger-cause Delta_Unemployment",
                 "Delta_Unemployment does not Granger-cause GDP_Growth"),
  F_stat  = round(c(gr_gdp$Granger$statistic, gr_unemp$Granger$statistic), 3),
  p_value = round(c(gr_gdp$Granger$p.value, gr_unemp$Granger$p.value), 4)
)
print(granger_tab, row.names = FALSE)
write.csv(granger_tab, "output/05_granger.csv", row.names = FALSE)
print(gr_gdp$Instant)

# ----------------------------------------------------------------------------
# 7. Orthogonalised impulse responses (Cholesky: GDP_Growth first)
# ----------------------------------------------------------------------------
set.seed(2026)
irf_gdp <- irf(var_model, impulse = "GDP_Growth", response = "Delta_Unemployment",
               n.ahead = IRF_HORIZON, ortho = TRUE, boot = TRUE,
               runs = IRF_RUNS, ci = 0.95, seed = 2026)
save_png("03_irf_gdp_to_unemployment.png",
         plot(irf_gdp, main = "Response of Unemployment to a GDP Growth Shock"))

set.seed(2026)
irf_unemp <- irf(var_model, impulse = "Delta_Unemployment", response = "GDP_Growth",
                 n.ahead = IRF_HORIZON, ortho = TRUE, boot = TRUE,
                 runs = IRF_RUNS, ci = 0.95, seed = 2026)
save_png("04_irf_unemployment_to_gdp.png",
         plot(irf_unemp, main = "Response of GDP Growth to an Unemployment Shock"))

# Point estimates and 95% bands for the headline IRF
irf_tab <- data.frame(
  horizon = 0:IRF_HORIZON,
  irf   = round(irf_gdp$irf$GDP_Growth[, 1], 3),
  lower = round(irf_gdp$Lower$GDP_Growth[, 1], 3),
  upper = round(irf_gdp$Upper$GDP_Growth[, 1], 3)
)
irf_tab$band_excludes_zero <- irf_tab$lower > 0 | irf_tab$upper < 0
print(irf_tab, row.names = FALSE)
write.csv(irf_tab, "output/06_irf_gdp_to_unemployment.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
# 8. Forecast error variance decomposition
# ----------------------------------------------------------------------------
fe <- fevd(var_model, n.ahead = IRF_HORIZON)
fevd_tab <- do.call(rbind, lapply(names(fe), function(v)
  data.frame(response = v, horizon = 1:IRF_HORIZON, round(100 * fe[[v]], 1))))
write.csv(fevd_tab, "output/06_fevd.csv", row.names = FALSE)
save_png("05_fevd.png", plot(fe), h = 1900)

# ----------------------------------------------------------------------------
# 9. Robustness
# ----------------------------------------------------------------------------

# 9a. Reverse Cholesky ordering (Delta_Unemployment first)
var_rev <- VAR(Y[, c("Delta_Unemployment", "GDP_Growth")], p = p, type = "const")
set.seed(2026)
irf_rev <- irf(var_rev, impulse = "GDP_Growth", response = "Delta_Unemployment",
               n.ahead = IRF_HORIZON, ortho = TRUE, boot = TRUE,
               runs = IRF_RUNS, ci = 0.95, seed = 2026)
save_png("06_irf_reverse_ordering.png",
         plot(irf_rev, main = "GDP Shock -> Unemployment: Reverse Ordering"))

# 9b. Pre-GFC subsample (1984Q1-2007Q4)
Y_preGFC   <- window(Y, end = c(2007, 4))
var_preGFC <- VAR(Y_preGFC, p = p, type = "const")
set.seed(2026)
irf_preGFC <- irf(var_preGFC, impulse = "GDP_Growth", response = "Delta_Unemployment",
                  n.ahead = IRF_HORIZON, ortho = TRUE, boot = TRUE,
                  runs = IRF_RUNS, ci = 0.95, seed = 2026)
save_png("07_irf_pre_gfc.png",
         plot(irf_preGFC, main = "GDP Shock -> Unemployment: Pre-GFC Sample"))

# Granger causality in the subsample (residuals pass normality here,
# so the standard F-test is used)
gr_pre <- data.frame(
  hypothesis = c("GDP_Growth does not Granger-cause Delta_Unemployment",
                 "Delta_Unemployment does not Granger-cause GDP_Growth"),
  F_stat  = round(c(causality(var_preGFC, cause = "GDP_Growth")$Granger$statistic,
                    causality(var_preGFC, cause = "Delta_Unemployment")$Granger$statistic), 3),
  p_value = round(c(causality(var_preGFC, cause = "GDP_Growth")$Granger$p.value,
                    causality(var_preGFC, cause = "Delta_Unemployment")$Granger$p.value), 4)
)
print(gr_pre, row.names = FALSE)
write.csv(gr_pre, "output/07_granger_pre_gfc.csv", row.names = FALSE)

# 9c. GFC impulse dummies (2008Q4, 2009Q1) on the full sample
tt  <- time(Y)
dum <- cbind(D08Q4 = as.numeric(abs(tt - 2008.75) < 0.01),
             D09Q1 = as.numeric(abs(tt - 2009.00) < 0.01))
stopifnot(sum(dum) == 2)
var_dum <- VAR(Y, p = p, type = "const", exogen = dum)

# 9d. Lag sensitivity of the headline IRF (VAR(1) to VAR(4))
lag_sens <- do.call(rbind, lapply(1:4, function(k) {
  m <- VAR(Y, p = k, type = "const")
  i <- irf(m, impulse = "GDP_Growth", response = "Delta_Unemployment",
           n.ahead = 6, ortho = TRUE, boot = FALSE)
  data.frame(lag = k, t(round(i$irf$GDP_Growth[1:6, 1], 3)))
}))
colnames(lag_sens) <- c("lag", paste0("h", 0:5))
print(lag_sens, row.names = FALSE)
write.csv(lag_sens, "output/07_lag_sensitivity.csv", row.names = FALSE)

# Diagnostics for all robustness models
models <- list(var_model, var_rev, var_preGFC, var_dum)
robustness_tab <- data.frame(
  specification = c("Baseline VAR(4)", "Reverse ordering",
                    "Pre-GFC subsample", "GFC dummies (full sample)"),
  bg_p_value = round(sapply(models, bg_pvalue), 4),
  jb_p_value = round(sapply(models, function(m)
    normality.test(m, multivariate.only = TRUE)$jb.mul$JB$p.value), 4),
  arch_p_value = round(sapply(models, function(m)
    arch.test(m, lags.multi = 2, multivariate.only = TRUE)$arch.mul$p.value), 4)
)
print(robustness_tab, row.names = FALSE)
write.csv(robustness_tab, "output/08_robustness_diagnostics.csv", row.names = FALSE)

# ----------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
cat("\nDone. Tables and figures are in the folder 'output'.\n")
