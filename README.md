# Output Growth and Unemployment Dynamics in the United States

### A bivariate VAR analysis of Okun's Law, 1984Q1–2019Q4

**Author:** Sameer Chawla  
**Programme:** MSc Economics, Gokhale Institute of Politics and Economics, Pune  
**Research type:** Applied time-series econometrics / reproducibility project  
**Sample:** US quarterly observations, 1984Q1–2019Q4

> **Main finding.** US real GDP growth and quarterly changes in unemployment exhibit a negative contemporaneous association. In an estimated bivariate VAR(3), lagged GDP growth helps predict unemployment-rate changes, whereas the reverse predictive hypothesis is not supported. The negative short-horizon orthogonalized impulse response is relatively stable across plausible lag orders, but its structural interpretation—and particularly the FEVD—is sensitive to recursive identification. Significant residual serial correlation and conditional heteroskedasticity are material limitations.

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

## 1. Research question and scope

**Question:** Do lagged changes in US real output growth contain information useful for predicting unemployment-rate changes, and how do these variables respond dynamically to identified innovations within a bivariate VAR?

This project evaluates the difference-form Okun association and then estimates a reduced-form VAR with both variables endogenous. The contributions are **empirical and pedagogical**, not claims of a novel estimator or a previously undocumented causal mechanism: the project contrasts a static benchmark with dynamic predictive tests, inspects residual adequacy, and measures sensitivity to lag specification, sample restriction and recursive identification.

A Granger-causality rejection is a **conditional predictive** result. It is not proof of an exogenous causal effect. Similarly, Cholesky-orthogonalized innovations are *identified by assumptions*, not uniquely recovered causal shocks.

## 2. Data and transformations

Source: [Federal Reserve Economic Data (FRED)](https://fred.stlouisfed.org/).

| Series | Source code | Native frequency | Treatment |
|---|---|---|---|
| Real gross domestic product | [`GDPC1`](https://fred.stlouisfed.org/series/GDPC1) | Quarterly, seasonally adjusted | 100 times the quarterly log difference |
| Civilian unemployment rate | [`UNRATE`](https://fred.stlouisfed.org/series/UNRATE) | Monthly, seasonally adjusted | Average the three monthly observations within each quarter; then first difference |

Let \(Q_t\) be quarterly real GDP and \(u_t\) the quarterly average unemployment rate (in percent). The variables are

\[
g_t = 100[\ln(Q_t)-\ln(Q_{t-1})], \qquad
\Delta u_t=u_t-u_{t-1}.
\]

**Units:** \(g_t\) is approximately percentage growth *per quarter*, **not** annualized; \(\Delta u_t\) is a change in the unemployment rate in **percentage points**. The first transformation requires 1983Q4 GDP, and the second requires 1983Q4 unemployment. The analytical sample has **144 observations** from 1984Q1 through 2019Q4; a full-sample VAR(3) uses **141 effective observations** after losing three lags.

The endpoint avoids including the extreme initial COVID-19 period; this is an explicit sample-scope choice, not a claim that the relation necessarily breaks in 2020. The beginning of the period roughly follows the start of the Great Moderation. Series from FRED are subject to data revisions: results may change slightly if re-downloaded at a later vintage.

| Variable | Mean | SD | Minimum | Maximum |
|---|---:|---:|---:|---:|
| GDP growth (%) | 0.6801 | 0.5723 | −2.2133 | 1.9359 |
| Change in unemployment (pp) | −0.0343 | 0.2691 | −0.6667 | 1.4000 |

The sample Pearson correlation is **−0.5814**; this does not establish a causal relationship.

## 3. Econometric specifications

### Static difference-form Okun benchmark

\[
\Delta u_t = \alpha + \beta g_t + e_t.
\]

Estimated by OLS, with **Newey–West HAC standard errors** (four lags; no prewhitening) reported for inference. This is a contemporaneous benchmark, not the VAR and not a causal structural equation.

### Reduced-form VAR

\[
\mathbf y_t = \mathbf c + \sum_{i=1}^{p} A_i\mathbf y_{t-i}+\boldsymbol\varepsilon_t,
\qquad
\mathbf y_t = \begin{bmatrix}g_t\\\Delta u_t\end{bmatrix},
\qquad E(\boldsymbol\varepsilon_t\boldsymbol\varepsilon_t')=\Sigma.
\]

The **version implemented in the uploaded R script uses VAR(3) with an intercept as its main specification**, with VAR(1) and VAR(4) as lag-sensitivity specifications. VAR coefficient matrices are estimated equation by equation using OLS. Dynamic stability is assessed using the eigenvalues of the companion matrix.

For orthogonalized impulse responses and FEVD, the baseline recursive ordering is **GDP growth first, change in unemployment second**. A Cholesky factor of \(\Sigma\) imposes a zero contemporaneous effect of the *second* identified innovation on GDP growth. This restriction does **not** follow from Granger causality and does not ensure the first innovation is a genuinely exogenous output shock.

## 4. Stationarity and lag selection

### Unit-root / stationarity evidence

ADF tests use an intercept and AIC-selected augmentation lags (maximum eight); PP tests use a constant with short and long bandwidths; KPSS tests assume level stationarity under the null and use both bandwidth choices.

| Test statistic | GDP growth | Δ unemployment | 5% critical value |
|---|---:|---:|---:|
| ADF | −5.1624 | −3.5448 | −2.8800 |
| PP, short bandwidth | −8.2749 | −5.3958 | −2.8817 |
| PP, long bandwidth | −8.6524 | −5.4091 | −2.8817 |
| KPSS, short bandwidth | **0.4861** | 0.1134 | 0.4630 |
| KPSS, long bandwidth | 0.3486 | 0.0827 | 0.4630 |

ADF and both PP versions reject their **unit-root nulls** at 5% for both transformed variables. Long-bandwidth KPSS does not reject its **stationarity null** for either series. **Short-bandwidth KPSS does reject level stationarity for GDP growth**; the discrepancy is reported rather than suppressed. Taken together, the tests support the working assumption of stationary transformed variables, with bandwidth sensitivity for GDP growth. They do **not** by themselves prove that either underlying level series has a particular integration order.

### Lag choice: a genuine specification trade-off

`VARselect` compares candidate lags through six on a common effective sample:

| Criterion | Selected lag |
|---|---:|
| AIC | 3 |
| FPE | 3 |
| HQ | 1 |
| SC/BIC | 1 |

The full-sample VAR(3) is used for the **AIC/FPE-preferred specification**, but it is **not fully adequate on residual diagnostics**: the 8-lag Breusch–Godfrey test rejects no serial correlation at 5% (p = 0.0443). The full-sample VAR(4) gives p = 0.1036 in the same nominal test, though that is not proof of zero residual dependence. The models have slightly different usable samples; separate common-sample comparisons align estimation to 138 quarters.

**Consequently VAR(3) should be described as the information-criterion baseline, not a uniquely validated final data-generating model.** VAR(4) is a substantive competing specification, not a cosmetic check.

## 5. Empirical findings

### 5.1 Static Okun regression

\[
\widehat{\Delta u_t}=0.1517-0.2734 g_t.
\]

| Statistic | Estimate |
|---|---:|
| GDP coefficient | −0.27339 |
| Newey–West HAC SE | 0.06879 |
| HAC t statistic | −3.9741 |
| HAC p value | 0.000112 |
| Approx. 95% HAC interval for GDP coefficient | [−0.4082, −0.1386] |
| OLS \(R^2\) | 0.3380 |

An extra percentage point of quarterly GDP growth is associated with a **0.273 percentage-point smaller quarterly change** in unemployment on average. It is not an estimate of an exogenous output intervention.

### 5.2 VAR(3) Granger predictability

Reported from the **full-sample VAR(3)** using `vars::causality(..., boot=TRUE, boot.runs=4999)` and the seeds shown in the code/log:

| Null hypothesis | F-statistic | Wild-bootstrap p value | Decision at 5% |
|---|---:|---:|---|
| Lagged GDP growth does not Granger-predict Δ unemployment | 8.0710 | **0.0224** | Reject |
| Lagged Δ unemployment does not Granger-predict GDP growth | 0.3852 | **0.7491** | Do not reject |

The evidence indicates asymmetric predictive content **conditional on the chosen variables, lag structure and sample**. It does not establish absence of economic feedback or structural causation. Classical F and HC1 joint Wald tests are additional checks in the script; HC1 covers heteroskedasticity under suitable assumptions but not arbitrary serial dependence.

### 5.3 VAR(3) impulse responses

For a positive **one-standard-deviation Cholesky-orthogonalized GDP-growth innovation**, the pointwise unemployment-change responses are:

| Horizon (quarters; 0 = impact) | Response (pp) | 95% bootstrap interval (pp) |
|---|---:|---|
| 0 | −0.0769 | [−0.1075, −0.0402] |
| 1 | −0.1004 | [−0.1346, −0.0597] |
| 2 | −0.0979 | [−0.1340, −0.0538] |
| 3 | −0.0455 | [−0.0806, −0.0023] |
| 4 | −0.0476 | [−0.0799, −0.0074] |
| 5 | −0.0318 | [−0.0613, 0.0044] |
| 12 | −0.0068 | [−0.0210, 0.0010] |

These intervals were reported from **1,999 standard residual-bootstrap runs**. They exclude zero at horizons 0–4 **pointwise**, and include zero from horizon 5 onward. They are *not simultaneous confidence bands*, and their nominal coverage may be unreliable under the detected conditional heteroskedasticity and residual dependence. Responses refer to **changes in unemployment**, not the unemployment level.

### 5.4 Forecast error variance decomposition (FEVD)

Under the GDP-first ordering, the contribution of GDP innovations to **forecast-error variance of Δ unemployment** is:

| Horizon | GDP innovation | Unemployment innovation |
|---|---:|---:|
| 1 quarter | 19.32% | 80.68% |
| 3 quarters | 44.33% | 55.67% |
| 12 quarters | **44.86%** | **55.14%** |

At twelve quarters, GDP-growth forecast error variance is approximately **98.63% own innovation** and **1.37% unemployment innovation** under the same ordering. These are *shares of uncertainty about forecasts*, **not** percentages of realized unemployment changes attributable to GDP.

## 6. Model diagnostics and limitations

The **full-sample VAR(3)** exhibits the following:

| Diagnostic | Null | Reported result | Interpretation |
|---|---|---|---|
| Breusch–Godfrey, 8 lags | No residual serial correlation | χ² = 46.791, p = **0.04427** | Rejected at 5% |
| Multivariate Jarque–Bera | Normal residuals | χ² = 75.751, p = **1.332 × 10⁻¹⁵** | Rejected |
| Multivariate ARCH-LM, 2 lags | No ARCH effects | χ² = 53.560, p = **2.149 × 10⁻⁵** | Rejected |
| Companion stability | All eigenvalue moduli below 1 | Maximum modulus **0.8075** | Dynamically stable |

The OLS-CUSUM processes in the included plots remain inside their displayed critical boundaries. This is a conditional stability diagnostic, not proof of no breaks of every kind. Notably, **dynamical stability is different from residual adequacy**.

The combined BG, JB and ARCH findings are important weaknesses, not technicalities. In particular, **the standard residual-bootstrap IRF intervals require caution**; a more demanding study could use inference explicitly robust to both conditional heteroskedasticity and remaining serial dependence, reassess lag length, and consider a richer information set. Coefficient consistency additionally requires appropriate conditional mean/exogeneity assumptions, rather than simply the absence of normality.

## 7. Robustness and identification sensitivity

**Lag specification.** A common-sample comparison aligns VAR(1), VAR(3) and VAR(4) to 138 usable quarters. Their initial GDP-to-Δ-unemployment IRFs are all negative and of broadly similar magnitude, although persistence and pointwise interval exclusions differ. Full-sample and common-sample diagnostic p-values must not be interchanged.

**Reverse Cholesky ordering.** Reordering variables leaves the reduced-form fit essentially the same after corresponding permutations, but changes the *definitions of orthogonal shocks*. For the full-sample VAR(3), the twelve-quarter GDP innovation contribution to Δ-unemployment FEVD changes from **44.86% (GDP first)** to **12.56% (unemployment first)**. Its impact effect on Δ unemployment becomes exactly zero in the reverse ordering by construction. Thus **the variance attribution is not identification-robust**, even though the GDP-growth-first IRF and reverse-ordered response share negative movements over some early horizons.

**Pre-global-financial-crisis sample.** The script re-estimates a VAR(3) on observations through 2007Q4 and plots its GDP-to-unemployment response. This checks sample sensitivity but does **not** prove the 2008–09 crisis is the sole source of all full-sample diagnostic rejections.

**Financial-crisis impulse indicators.** The script adds separate 2008Q4 and 2009Q1 impulse dummies as a descriptive intervention sensitivity check. These do not by themselves provide an economically identified crisis model or guarantee removal of residual misspecification.

**Omitted-variable and policy confounding.** Labour-force participation, productivity, structural labour-market changes and macroeconomic policy are not explicitly included. A two-variable Cholesky system cannot isolate common third-factor shocks or prove a policy-relevant exogenous output innovation.

**Data vintage.** `GDPC1` is regularly revised. Re-running FRED downloads at a later date can alter numerical results even with unchanged code and seeds. A fixed-vintage dataset would be needed for strict numerical replication independent of live FRED revisions.

## 8. Figures

The following paths are those generated by the **uploaded VAR(3) R script**. They display only once the script has created the `output/` directory and the PNGs have been committed to GitHub. Do not substitute VAR(4) pictures under VAR(3) captions.

**Figure 1.** Quarterly real GDP growth and change in unemployment, 1984Q1–2019Q4.

![Quarterly GDP growth and unemployment-rate changes](output/01_series.png)

**Figure 2.** Static difference-form Okun association (OLS line; HAC inference reported above).

![Static Okun scatterplot](output/08_okun_scatter.png)

**Figure 3.** VAR(3) orthogonalized GDP-growth innovation → change in unemployment, pointwise 95% bootstrap intervals.

![GDP-to-unemployment VAR(3) IRF](output/03_irf_gdp_to_unemployment.png)

**Figure 4.** Baseline VAR(3) forecast-error variance decomposition.

![VAR(3) forecast error variance decomposition](output/05_fevd.png)

**Figure 5.** Common-sample lag-order sensitivity of the GDP-to-unemployment IRF.

![Lag sensitivity](output/09_lag_sensitivity_irf.png)

<details>
<summary>Additional figures: identification, sample sensitivity and diagnostics</summary>

**Figure A1.** OLS-CUSUM diagnostics for both VAR equations.

![CUSUM](output/02_cusum.png)

**Figure A2.** GDP response to an identified unemployment-change innovation (GDP ordered first).

![Reverse-direction IRF](output/04_irf_unemployment_to_gdp.png)

**Figure A3.** GDP innovation → unemployment change, reversed Cholesky ordering.

![Reverse Cholesky ordering](output/06_irf_reverse_ordering.png)

**Figure A4.** GDP innovation → unemployment change, pre-GFC subsample.

![Pre-GFC subsample](output/07_irf_pre_gfc.png)

</details>

> **Presentation note:** The current package-default FEVD figures can place the legend over the bars. This is a plotting/layout issue, not an econometric result; improve the exported FEVD layout before using it in a formal PDF submission. IRF confidence intervals are pointwise, and the bootstrap seed/replication count should be stated in any independent presentation.

## 9. Reproducibility

### Repository structure

```text
US-OKUNS-LAW-VAR/
├── README.md
├── Okun_Law_Bivariate_VAR (1).R   # currently uploaded script name; rename if desired
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
    └── table_*.csv
```

### Software and execution

- Tested by the author using **R 4.5.0 on Windows** for the constituent analyses; the exact consolidated repository script should also be run from a fresh R session as a final end-to-end check.
- R packages used by the script: `quantmod`, `zoo`, `urca`, `vars`, `sandwich`, `lmtest` (and their dependencies).
- Internet access required for the live FRED downloads.
- Bootstrap settings in the currently uploaded VAR(3) script: **4,999 Granger runs**, **1,999 IRF runs**; changing these changes Monte Carlo precision and possibly the reported p-values/interval endpoints.

```r
install.packages(c("quantmod", "zoo", "urca", "vars", "sandwich", "lmtest"))
source("Okun_Law_Bivariate_VAR (1).R")
```

Run from the repository root. The script creates `output/`, writes the CSV tables and PNG figures, and prints estimates and tests to the R console. If you rename the script to `Okun_Law_Bivariate_VAR.R`, change the `source()` command accordingly. If pictures do not render on GitHub, confirm the `output/` directory has actually been committed (not only generated locally).

**Reproducibility status:** Results in this README were transcribed from supplied R console runs and accompanying plots; the exact current consolidated script has **not been independently executed end to end in this audit**. The repository should not claim byte-for-byte reproducibility until a clean run reproduces its numerical tables and all referenced figures. For stringent archival replication, record a fixed FRED/ALFRED vintage and package versions.

## 10. References

- Dickey, D. A., & Fuller, W. A. (1979). Distribution of the estimators for autoregressive time series with a unit root. *Journal of the American Statistical Association*, 74(366), 427–431.
- Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica*, 37(3), 424–438.
- Kwiatkowski, D., Phillips, P. C. B., Schmidt, P., & Shin, Y. (1992). Testing the null hypothesis of stationarity against the alternative of a unit root. *Journal of Econometrics*, 54(1–3), 159–178.
- Newey, W. K., & West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. *Econometrica*, 55(3), 703–708.
- Okun, A. M. (1962). Potential GNP: Its measurement and significance. *Proceedings of the Business and Economic Statistics Section*, American Statistical Association.
- Phillips, P. C. B., & Perron, P. (1988). Testing for a unit root in time series regression. *Biometrika*, 75(2), 335–346.
- Pfaff, B. (2008). VAR, SVAR and SVEC models: Implementation within R package `vars`. *Journal of Statistical Software*, 27(4), 1–32.
- U.S. Bureau of Economic Analysis. Real Gross Domestic Product [`GDPC1`]. Federal Reserve Bank of St. Louis, FRED. https://fred.stlouisfed.org/series/GDPC1
- U.S. Bureau of Labor Statistics. Unemployment Rate [`UNRATE`]. Federal Reserve Bank of St. Louis, FRED. https://fred.stlouisfed.org/series/UNRATE

---

*Research code and empirical results are provided for academic examination and replication. The analysis concerns predictive relationships and conditional identified innovations; it should not be interpreted as establishing structural causal effects.*
