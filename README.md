# Output Growth and Unemployment Dynamics in the US

### A bivariate VAR analysis of Okun's Law, 1984Q1–2019Q4

![R](https://img.shields.io/badge/R-4.5.0-276DC3?logo=r&logoColor=white)
![Data](https://img.shields.io/badge/data-FRED-1f6feb)
![Sample](https://img.shields.io/badge/sample-1984Q1--2019Q4-informational)
![Model](https://img.shields.io/badge/model-VAR(3)-success)

**Author:** Sameer Chawla · MSc Economics, Gokhale Institute of Politics and Economics, Pune
**Type:** Applied time-series econometrics · reproducible R project

<p align="center">
  <img src="output/03_irf_gdp_to_unemployment.png" width="48%" alt="Impulse response of unemployment change to a GDP growth shock">
  <img src="output/08_okun_scatter.png" width="48%" alt="Difference-form Okun scatterplot with OLS line">
</p>

> **Main finding.** Faster US GDP growth goes with a *smaller* change in the unemployment rate, and lagged GDP growth helps predict unemployment changes while the reverse does not. A one-standard-deviation GDP shock lowers the unemployment-rate change by about **0.08 pp on impact and 0.10 pp after one quarter**, fading over roughly a year. The shape of the response is stable across lag orders, but the **variance decomposition depends on the Cholesky ordering**, and residual diagnostics reveal real limitations (see [Section 6](#6-model-diagnostics-and-limitations)).

---

## Results at a glance

| Question | Result |
|---|---|
| Static Okun relation | $\widehat{\Delta u_t} = 0.1517 - 0.2734\,g_t$; HAC $t = -3.97$; $R^2 = 0.338$ |
| Does GDP growth Granger-predict Δ unemployment? | **Yes**: F = 8.07, bootstrap p = 0.022 |
| Does Δ unemployment Granger-predict GDP growth? | **No**: F = 0.39, bootstrap p = 0.749 |
| Peak IRF response (GDP shock → Δ unemployment) | −0.100 pp at horizon 1 (95% CI −0.135 to −0.060) |
| 12-quarter FEVD share of GDP shocks in Δ unemployment | 44.9% (GDP first) vs **12.6%** (unemployment first) |
| Model stability | Stable (largest companion eigenvalue modulus 0.808) |
| Residual diagnostics | Serial correlation, non-normality and ARCH effects detected |

## Contents

1. [Research question and scope](#1-research-question-and-scope)
2. [Data and transformations](#2-data-and-transformations)
3. [Econometric specifications](#3-econometric-specifications)
4. [Stationarity and lag selection](#4-stationarity-and-lag-selection)
5. [Empirical findings](#5-empirical-findings)
6. [Model diagnostics and limitations](#6-model-diagnostics-and-limitations)
7. [Robustness and identification sensitivity](#7-robustness-and-identification-sensitivity)
8. [Figures](#8-figures)
9. [Reproducibility](#9-reproducibility)
10. [References](#10-references)

---

## 1. Research question and scope

**Question:** Do lagged changes in US real output growth help predict unemployment-rate changes, and how do the two variables respond dynamically to identified innovations in a bivariate VAR?

The project estimates the difference-form Okun relationship, then a reduced-form VAR with both variables endogenous. The contribution is **empirical and pedagogical**: it contrasts a static benchmark with dynamic predictive tests, inspects residual adequacy, and measures sensitivity to lag length, sample and recursive identification. It does not claim a new estimator or a newly discovered causal mechanism.

Two interpretive cautions apply throughout:

- A Granger-causality rejection is a **conditional predictive** result, not proof of an exogenous causal effect.
- Cholesky-orthogonalized innovations are **identified by assumption**, not uniquely recovered structural shocks.

## 2. Data and transformations

Source: [Federal Reserve Economic Data (FRED)](https://fred.stlouisfed.org/).

| Series | FRED code | Native frequency | Treatment |
|---|---|---|---|
| Real GDP | [`GDPC1`](https://fred.stlouisfed.org/series/GDPC1) | Quarterly, seasonally adjusted | 100 × quarterly log difference |
| Civilian unemployment rate | [`UNRATE`](https://fred.stlouisfed.org/series/UNRATE) | Monthly, seasonally adjusted | Average the three monthly values per quarter, then first-difference |

With $Q_t$ as quarterly real GDP and $u_t$ as the quarterly average unemployment rate (percent):

$$
g_t = 100\left[\ln Q_t - \ln Q_{t-1}\right], \qquad \Delta u_t = u_t - u_{t-1}.
$$

**Units.** $g_t$ is approximately percentage growth *per quarter* (not annualized). $\Delta u_t$ is in **percentage points**.

**Sample.** 144 observations, 1984Q1–2019Q4 (the transformations need 1983Q4 values). A VAR(3) uses **141 effective observations** after losing three lags. The sample starts near the beginning of the Great Moderation and ends before COVID-19; the endpoint is a scope choice, not a claim that the relationship breaks in 2020. FRED data are revised, so re-downloading later may change results slightly.

| Variable | Mean | SD | Min | Max |
|---|---:|---:|---:|---:|
| GDP growth (%) | 0.6801 | 0.5723 | −2.2133 | 1.9359 |
| Change in unemployment (pp) | −0.0343 | 0.2691 | −0.6667 | 1.4000 |

The sample correlation is **−0.5814**, which describes association and does not establish causation.

![Quarterly GDP growth and unemployment-rate changes](output/01_series.png)

## 3. Econometric specifications

### Static difference-form Okun benchmark

$$
\Delta u_t = \alpha + \beta\, g_t + e_t
$$

Estimated by OLS with **Newey–West HAC standard errors** (4 lags, no prewhitening). This is a contemporaneous benchmark, not a structural equation.

### Reduced-form VAR

$$
\mathbf{y}_t = \mathbf{c} + \sum_{i=1}^{p} A_i \mathbf{y}_{t-i} + \boldsymbol{\varepsilon}_t,
\qquad
\mathbf{y}_t = \begin{bmatrix} g_t \\ \Delta u_t \end{bmatrix},
\qquad
E(\boldsymbol{\varepsilon}_t \boldsymbol{\varepsilon}_t') = \Sigma
$$

The main specification is **VAR(3) with an intercept**, estimated equation by equation with OLS. VAR(1) and VAR(4) serve as lag-sensitivity comparisons. Stability is checked with the eigenvalues of the companion matrix.

**Identification.** For orthogonalized IRFs and FEVD the baseline ordering is **GDP growth first, Δ unemployment second**. The Cholesky factor of $\Sigma$ then imposes that the second innovation has no contemporaneous effect on GDP growth. This restriction does not follow from Granger causality and does not make the first innovation a genuinely exogenous output shock.

## 4. Stationarity and lag selection

### Unit-root and stationarity tests

ADF uses an intercept with AIC-selected lags (max 8). PP uses a constant with short and long bandwidths. KPSS takes level stationarity as the null and uses both bandwidths.

| Test | GDP growth | Δ unemployment | 5% critical value |
|---|---:|---:|---:|
| ADF | −5.1624 | −3.5448 | −2.8800 |
| PP, short bandwidth | −8.2749 | −5.3958 | −2.8817 |
| PP, long bandwidth | −8.6524 | −5.4091 | −2.8817 |
| KPSS, short bandwidth | **0.4861** | 0.1134 | 0.4630 |
| KPSS, long bandwidth | 0.3486 | 0.0827 | 0.4630 |

ADF and both PP versions reject a unit root for both variables. Long-bandwidth KPSS does not reject stationarity for either. **Short-bandwidth KPSS does reject stationarity for GDP growth**, a discrepancy reported rather than hidden. Overall the tests support treating both transformed series as stationary, with some bandwidth sensitivity for GDP growth.

### Lag selection is a genuine trade-off

| Criterion | Selected lag |
|---|---:|
| AIC | 3 |
| FPE | 3 |
| HQ | 1 |
| SC/BIC | 1 |

VAR(3) follows AIC/FPE but is **not fully adequate on residual diagnostics**: the 8-lag Breusch–Godfrey test rejects no serial correlation at 5% (p = 0.0443). The full-sample VAR(4) gives p = 0.1036, which is not proof of zero residual dependence. Because models have slightly different usable samples, common-sample comparisons align estimation to **138 quarters**.

> VAR(3) is best described as the **information-criterion baseline**, not a uniquely validated data-generating model. VAR(4) is a substantive competitor.

## 5. Empirical findings

### 5.1 Static Okun regression

$$
\widehat{\Delta u_t} = 0.1517 - 0.2734\, g_t
$$

| Statistic | Estimate |
|---|---:|
| GDP coefficient | −0.27339 |
| Newey–West HAC SE | 0.06879 |
| HAC t statistic | −3.9741 |
| HAC p value | 0.000112 |
| Approx. 95% HAC interval | [−0.4082, −0.1386] |
| OLS $R^2$ | 0.3380 |

One extra percentage point of quarterly GDP growth is associated with a **0.273 pp smaller** quarterly change in unemployment, on average. This is an association, not the effect of an exogenous output intervention.

![Static Okun scatterplot](output/08_okun_scatter.png)

### 5.2 VAR(3) Granger predictability

Computed on the full-sample VAR(3) with `vars::causality(..., boot = TRUE, boot.runs = 4999)`.

| Null hypothesis | F-statistic | Wild-bootstrap p | Decision (5%) |
|---|---:|---:|---|
| Lagged GDP growth does not Granger-predict Δ unemployment | 8.0710 | **0.0224** | Reject |
| Lagged Δ unemployment does not Granger-predict GDP growth | 0.3852 | **0.7491** | Do not reject |

The evidence points to asymmetric predictive content **given these variables, this lag structure and this sample**. It does not rule out economic feedback or establish structural causation. Classical F and HC1 joint Wald tests are included in the script as additional checks; HC1 addresses heteroskedasticity but not arbitrary serial dependence.

### 5.3 VAR(3) impulse responses

Response of Δ unemployment to a positive **one-standard-deviation Cholesky-orthogonalized GDP-growth innovation**:

| Horizon (quarters; 0 = impact) | Response (pp) | 95% bootstrap interval (pp) |
|---|---:|---|
| 0 | −0.0769 | [−0.1075, −0.0402] |
| 1 | −0.1004 | [−0.1346, −0.0597] |
| 2 | −0.0979 | [−0.1340, −0.0538] |
| 3 | −0.0455 | [−0.0806, −0.0023] |
| 4 | −0.0476 | [−0.0799, −0.0074] |
| 5 | −0.0318 | [−0.0613, 0.0044] |
| 12 | −0.0068 | [−0.0210, 0.0010] |

Intervals come from **1,999 residual-bootstrap runs**. They exclude zero **pointwise** at horizons 0–4 and include zero from horizon 5 onward. They are not simultaneous bands, and their coverage may be unreliable given the ARCH effects and residual dependence found in Section 6. Responses are in **changes** in unemployment, not its level.

![VAR(3) impulse response](output/03_irf_gdp_to_unemployment.png)

### 5.4 Forecast error variance decomposition (GDP-first ordering)

Share of **Δ unemployment** forecast-error variance:

| Horizon | GDP innovation | Unemployment innovation |
|---|---:|---:|
| 1 quarter | 19.32% | 80.68% |
| 3 quarters | 44.33% | 55.67% |
| 12 quarters | **44.86%** | **55.14%** |

For GDP growth itself, about **98.63%** of the 12-quarter forecast-error variance is its own innovation and **1.37%** comes from unemployment. These are shares of *forecast uncertainty*, not the percentage of realized unemployment changes caused by GDP.

## 6. Model diagnostics and limitations

Full-sample VAR(3):

| Diagnostic | Null hypothesis | Result | Reading |
|---|---|---|---|
| Breusch–Godfrey (8 lags) | No serial correlation | χ² = 46.791, p = **0.0443** | Rejected at 5% |
| Multivariate Jarque–Bera | Normal residuals | χ² = 75.751, p = **1.33 × 10⁻¹⁵** | Rejected |
| Multivariate ARCH-LM (2 lags) | No ARCH effects | χ² = 53.560, p = **2.15 × 10⁻⁵** | Rejected |
| Companion stability | All moduli < 1 | Max modulus **0.8075** | Stable |

OLS-CUSUM paths stay inside their critical bounds (Figure A1). That is a conditional stability check, not proof of no breaks of any kind, and **dynamic stability is not the same as residual adequacy**.

The serial correlation, non-normality and ARCH results are substantive weaknesses. In particular, **the residual-bootstrap IRF intervals should be read with caution**. A more demanding study could use inference robust to conditional heteroskedasticity and remaining serial dependence, revisit lag length, and widen the information set.

## 7. Robustness and identification sensitivity

**Lag order.** VAR(1), VAR(3) and VAR(4) are aligned to a common 138-quarter sample. All three give negative initial responses of broadly similar size; persistence and pointwise interval exclusions differ. Full-sample and common-sample diagnostic p-values must not be mixed.

![Lag-order sensitivity](output/09_lag_sensitivity_irf.png)

**Reverse Cholesky ordering.** Reordering leaves the reduced-form fit unchanged but redefines the orthogonal shocks. The 12-quarter GDP contribution to Δ unemployment FEVD moves from **44.86% (GDP first)** to **12.56% (unemployment first)**, and the impact response becomes exactly zero by construction. **Variance attribution is therefore not robust to identification.**

**Pre-GFC sample.** A VAR(3) re-estimated on data through 2007Q4 checks sample sensitivity. It does not show the 2008–09 crisis is the only source of the full-sample diagnostic failures.

**Financial-crisis dummies.** Separate impulse dummies for 2008Q4 and 2009Q1 provide a descriptive intervention check, not an economically identified crisis model.

**Omitted variables.** Labour-force participation, productivity, structural labour-market change and policy are not modelled. A two-variable system cannot separate common third-factor shocks.

**Data vintage.** `GDPC1` is revised regularly; strict numerical replication needs a fixed vintage (e.g. from ALFRED).

## 8. Figures

| # | File | Description |
|---|---|---|
| 1 | `01_series.png` | GDP growth and Δ unemployment, 1984Q1–2019Q4 |
| 2 | `08_okun_scatter.png` | Static Okun scatter with OLS line |
| 3 | `03_irf_gdp_to_unemployment.png` | VAR(3) IRF, GDP shock → Δ unemployment |
| 4 | `05_fevd.png` | VAR(3) forecast-error variance decomposition |
| 5 | `09_lag_sensitivity_irf.png` | IRF across VAR(1), VAR(3), VAR(4) on a common sample |

**Figure 4. Variance decomposition (GDP first).**

![VAR(3) FEVD](output/05_fevd.png)

<details>
<summary><b>Additional figures: diagnostics, identification and sample sensitivity</b></summary>

**Figure A1. OLS-CUSUM stability plots for both equations.**

![CUSUM](output/02_cusum.png)

**Figure A2. VAR residuals (GDP growth vs Δ unemployment).**

![VAR residuals](output/10_var_residuals.png)

**Figure A3. Response of GDP growth to a Δ-unemployment innovation (GDP ordered first).**

![Reverse-direction IRF](output/04_irf_unemployment_to_gdp.png)

**Figure A4. GDP innovation → Δ unemployment under the reversed ordering (impact response is zero by construction).**

![Reverse ordering](output/06_irf_reverse_ordering.png)

**Figure A5. GDP innovation → Δ unemployment, pre-GFC subsample (through 2007Q4).**

![Pre-GFC subsample](output/07_irf_pre_gfc.png)

</details>

## 9. Reproducibility

### Repository structure

```text
US-OKUNS-LAW-VAR/
├── README.md
├── Okun_Law_Bivariate_VAR.R
├── LICENSE
└── output/
    ├── 01_series.png
    ├── 02_cusum.png
    ├── 03_irf_gdp_to_unemployment.png
    ├── 04_irf_unemployment_to_gdp.png
    ├── 05_fevd.png
    ├── 06_irf_reverse_ordering.png
    ├── 07_irf_pre_gfc.png
    ├── 08_okun_scatter.png
    ├── 09_lag_sensitivity_irf.png
    ├── 10_var_residuals.png
    └── table_*.csv
```

### Quick start

```r
install.packages(c("quantmod", "zoo", "urca", "vars", "sandwich", "lmtest"))
source("Okun_Law_Bivariate_VAR.R")   # run from the repository root
```

The script downloads data from FRED (internet required), creates `output/`, writes the CSV tables and PNG figures, and prints estimates and tests to the console.

| Item | Detail |
|---|---|
| R version | 4.5.0 (Windows) |
| Packages | `quantmod`, `zoo`, `urca`, `vars`, `sandwich`, `lmtest` |
| Granger bootstrap | 4,999 runs |
| IRF bootstrap | 1,999 runs (main VAR(3) IRF) |
| Data | Live FRED download (not a fixed vintage) |

Changing bootstrap settings changes Monte Carlo precision and may shift reported p-values and interval endpoints.

**Reproducibility status.** Reported numbers come from R console runs of the VAR(3) analysis. Before treating the repository as an archival replication package, run the consolidated script once in a fresh R session, confirm it regenerates every table and figure above, and record package versions (`sessionInfo()`) and a fixed FRED/ALFRED vintage.

## 10. References

- Dickey, D. A., & Fuller, W. A. (1979). Distribution of the estimators for autoregressive time series with a unit root. *Journal of the American Statistical Association*, 74(366), 427–431.
- Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*, 37(3), 424–438.
- Kwiatkowski, D., Phillips, P. C. B., Schmidt, P., & Shin, Y. (1992). Testing the null hypothesis of stationarity against the alternative of a unit root. *Journal of Econometrics*, 54(1–3), 159–178.
- Newey, W. K., & West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. *Econometrica*, 55(3), 703–708.
- Okun, A. M. (1962). Potential GNP: Its measurement and significance. *Proceedings of the Business and Economic Statistics Section*, American Statistical Association.
- Pfaff, B. (2008). VAR, SVAR and SVEC models: Implementation within R package `vars`. *Journal of Statistical Software*, 27(4), 1–32.
- Phillips, P. C. B., & Perron, P. (1988). Testing for a unit root in time series regression. *Biometrika*, 75(2), 335–346.
- U.S. Bureau of Economic Analysis. Real Gross Domestic Product [`GDPC1`]. FRED. https://fred.stlouisfed.org/series/GDPC1
- U.S. Bureau of Labor Statistics. Unemployment Rate [`UNRATE`]. FRED. https://fred.stlouisfed.org/series/UNRATE

---

*Research code and results are provided for academic examination and replication. The analysis concerns predictive relationships and assumption-identified innovations and should not be read as establishing structural causal effects.*
