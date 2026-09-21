ynamic Okun's Law in the United States

A Bivariate VAR Analysis of Real GDP Growth and Changes in Unemployment, 1984Q1–2019Q4






This project studies the short-run dynamic relationship between U.S. real GDP growth and changes in the unemployment rate using a bivariate Vector Autoregression (VAR). The goal is not to claim structural causality, but to test whether output growth contains predictive information for subsequent unemployment changes and to examine how the two variables respond dynamically to innovations in the system.

The project is designed as an MSc-level applied econometrics exercise with emphasis on transparent model specification, stationarity testing, lag selection, residual diagnostics, Granger causality, impulse responses, forecast-error variance decomposition, and robustness checks.

Contents

Research question

Data and variables

Econometric strategy

Main results

Impulse-response analysis

Variance decomposition

Diagnostics and robustness

Limitations

Repository structure

How to run

References

Research question

Do changes in U.S. real economic activity systematically precede changes in unemployment, and how persistent is that relationship?

The project focuses on a dynamic version of Okun's Law. A static Okun regression summarizes contemporaneous co-movement between output and unemployment. A VAR allows both variables to depend on their own lags and on past values of the other variable.

The estimated reduced-form system is:

Y_t = c + A_1 Y_(t-1) + ... + A_p Y_(t-p) + e_t

where

Y_t = [ GDP_Growth_t , Delta_Unemployment_t ]'

The two variables are treated as jointly endogenous.

For orthogonalized impulse responses, the baseline recursive ordering is:

GDP_Growth  ->  Delta_Unemployment

This means GDP-growth innovations may affect unemployment within the same quarter under the baseline Cholesky identification. Because that is an identifying assumption rather than a fact established by the data, the reverse ordering is also reported as a robustness check.

Data and variables

Source: Federal Reserve Economic Data (FRED), downloaded directly in R using quantmod.

Series

FRED code

Frequency in source

Use in project

Real Gross Domestic Product

GDPC1

Quarterly

Converted to quarter-on-quarter log growth

Civilian Unemployment Rate

UNRATE

Monthly

Aggregated to quarterly mean, then first-differenced

Sample: 1984Q1–2019Q4
Observations: 144 quarters

The sample ends in 2019Q4 so that the unusually large COVID-19 shock does not dominate the covariance structure of a linear VAR.

Variable construction

GDP_Growth_t = 100 x [ ln(GDP_t) - ln(GDP_(t-1)) ]

Delta_Unemployment_t = Unemployment_t - Unemployment_(t-1)

GDP_Growth is therefore quarter-on-quarter real GDP log growth in percent.
Delta_Unemployment is the quarterly change in the unemployment rate in percentage points.

Descriptive statistics

Variable

Mean

SD

Minimum

Maximum

GDP_Growth

0.680

0.572

-2.213

1.936

Delta_Unemployment

-0.034

0.269

-0.667

1.400

Time-series plot



The large 2008–09 movements are visible in both series and later matter for the residual normality and heteroskedasticity diagnostics.

Econometric strategy

1. Stationarity

The variables entering the VAR are tested directly.

Test

GDP_Growth

Delta_Unemployment

Interpretation

ADF statistic

-5.206

-3.906

Reject unit-root null at 5%

ADF 5% critical value

-2.88

-2.88



KPSS, long-lag statistic

0.349

0.083

Fail to reject stationarity at 5%

KPSS 5% critical value

0.463

0.463



The ADF and long-lag KPSS results therefore support treating both transformed variables as stationary, I(0).

A short-lag KPSS specification for GDP growth is borderline, so stationarity is treated as strongly supported by the combined evidence rather than as something mechanically "proved" by one test.

2. Lag selection

Information criteria do not agree:

Criterion

Preferred lag

Schwarz / BIC

1

Hannan-Quinn

1

AIC

3

FPE

3

The parsimonious VAR(1) leaves residual autocorrelation:

Breusch-Godfrey p-value for VAR(1) = 0.013

The lag order is therefore increased sequentially. VAR(4) is the shortest specification for which the 8-lag Breusch-Godfrey test fails to reject residual serial correlation:

Final lag order = 4
Breusch-Godfrey p-value = 0.104

This is best viewed as diagnostic respecification rather than a claim that BIC itself selected four lags.

3. Stability

The largest companion-matrix root modulus is:

0.706

All roots lie inside the unit circle, so the final VAR(4) is dynamically stable.

Main results

Granger-causality tests

Bootstrap inference is used for the Granger tests.

Null hypothesis

F statistic

Bootstrap p-value

Result

GDP growth does not Granger-cause change in unemployment

6.164

0.023

Reject

Change in unemployment does not Granger-cause GDP growth

0.282

0.863

Do not reject

There is also strong contemporaneous dependence between the reduced-form innovations:

Instantaneous-causality test: p < 0.001
Residual innovation correlation: approximately -0.45

Interpretation

The evidence is asymmetric in predictive terms:

Past GDP growth contains statistically significant information for forecasting subsequent changes in unemployment, conditional on the VAR's own lag structure. The reverse predictive relationship is not supported.

This is a Granger-predictive result, not proof that GDP growth is an exogenous structural cause of unemployment.

Impulse-response analysis

Baseline ordering: GDP growth first

A positive orthogonalized GDP-growth innovation produces a negative response in the change in unemployment. The response is strongest in the first few quarters and then gradually decays toward zero.



Reverse experiment

The response of GDP growth to an unemployment innovation is negative in the orthogonalized system, although the Granger test does not support a lagged predictive effect running from unemployment to GDP growth.



Identification robustness

Reversing the Cholesky ordering preserves the negative short-run sign of the GDP-growth-shock response of unemployment, but the magnitude changes.



This is an important qualification: the qualitative short-run relationship is reasonably robust, while exact structural magnitudes remain identification-dependent.

Variance decomposition

Under the baseline Cholesky ordering:

Response

Horizon

Own innovation

Other innovation

GDP_Growth

1

100.0%

0.0%

GDP_Growth

12

98.7%

1.3%

Delta_Unemployment

1

79.5%

20.5%

Delta_Unemployment

12

55.8%

44.2%



The baseline FEVD suggests that GDP-growth innovations account for a substantial share of unemployment forecast-error variance at longer horizons.

However, this number is sensitive to recursive ordering. Under the reverse ordering, the GDP-growth share of unemployment forecast-error variance at horizon 12 falls substantially. For that reason, the FEVD is interpreted as an identification-dependent decomposition rather than as a structural fact.

Diagnostics and robustness

Baseline VAR(4)

Diagnostic

Result

Interpretation

Breusch-Godfrey LM

p = 0.104

No evidence of remaining serial correlation at the tested horizon

Largest VAR root

0.706

Dynamically stable

Jarque-Bera

p < 0.001

Residual normality rejected

Multivariate ARCH-LM

p < 0.001

Conditional homoskedasticity rejected

OLS-CUSUM

Within critical bands

No formal evidence of coefficient instability

CUSUM stability



The CUSUM processes remain within the critical boundaries. This does not eliminate all concerns about crisis-period instability, but it provides no formal CUSUM evidence of parameter instability.

Pre-GFC robustness

The model is re-estimated over 1984Q1–2007Q4.



The negative GDP-growth-shock response of unemployment is preserved. The shorter pre-GFC specification also has substantially cleaner residual diagnostics, suggesting that much of the full-sample non-normality and conditional heteroskedasticity is associated with the Global Financial Crisis period.

Crisis-dummy and lag-sensitivity checks

The full script additionally evaluates:

impulse dummies around the Global Financial Crisis;

alternative VAR lag orders;

reverse recursive ordering;

the pre-GFC subsample.

These checks are used to assess whether the headline result is being driven mechanically by one specific lag choice, one crisis episode, or one Cholesky ordering.

Limitations

The project is intentionally small and transparent, but several limitations matter.

Residual heteroskedasticity and non-normality.
The full-sample VAR rejects multivariate normality and homoskedasticity. OLS point estimates remain useful under appropriate conditions, but conventional Gaussian inference and standard residual-bootstrap IRF bands should be interpreted cautiously.

Recursive identification.
Orthogonalized IRFs and FEVD require a Cholesky ordering. Reversing the ordering preserves the qualitative short-run IRF sign but changes the magnitude and materially changes the FEVD. Structural interpretations are therefore deliberately limited.

Bivariate information set.
The model contains only output growth and unemployment changes. Productivity, labour-force participation, monetary policy, fiscal policy and other common drivers are omitted. The model is therefore best viewed as a compact reduced-form study of dynamic Okun-type co-movement.

Sample choice.
The sample excludes 2020 onward. This avoids allowing COVID-19 to dominate estimation but means the model is not intended to describe the pandemic or post-pandemic labour market.

Granger causality is not structural causality.
A statistically significant Granger test means one variable contains incremental predictive information for another conditional on the model. It does not by itself establish an exogenous causal mechanism.

Repository structure

US-OKUNS-LAW-VAR/
|
|-- README.md
|-- Okun_Law_Bivariate_VAR.R
|-- LICENSE
|
`-- output/
    |-- 01_series.png
    |-- 02_cusum.png
    |-- 03_irf_gdp_to_unemployment.png
    |-- 04_irf_unemployment_to_gdp.png
    |-- 05_fevd.png
    |-- 06_irf_reverse_ordering.png
    `-- 07_irf_pre_gfc.png

Important: the image files must actually exist inside the repository at these exact paths. GitHub paths are case-sensitive.

How to run

Requirements

R 4.5 or later

Internet connection for FRED downloads

Packages:

quantmod

zoo

urca

vars

Run

Clone or download the repository, set the repository as the working directory, and run:

source("Okun_Law_Bivariate_VAR.R")

The script downloads the FRED data, constructs the quarterly variables, performs the stationarity and diagnostic tests, estimates the VAR, produces the causality/IRF/FEVD analysis, and runs the robustness checks.

Reproducibility notes

FRED macroeconomic series can be revised over time. Exact numerical results may therefore change slightly when the script is rerun with later data vintages.

Random seeds are set before bootstrap procedures to improve reproducibility.

References

Breusch, T. S. (1978). Testing for autocorrelation in dynamic linear models. Australian Economic Papers, 17(31), 334–355.

Dickey, D. A., & Fuller, W. A. (1979). Distribution of the estimators for autoregressive time series with a unit root. Journal of the American Statistical Association, 74(366a), 427–431.

Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. Econometrica, 37(3), 424–438.

Kwiatkowski, D., Phillips, P. C. B., Schmidt, P., & Shin, Y. (1992). Testing the null hypothesis of stationarity against the alternative of a unit root. Journal of Econometrics, 54(1–3), 159–178.

Lütkepohl, H. (2005). New Introduction to Multiple Time Series Analysis. Springer.

Okun, A. M. (1962). Potential GNP: Its Measurement and Significance. Proceedings of the Business and Economic Statistics Section, American Statistical Association, 98–104.

Sims, C. A. (1980). Macroeconomics and reality. Econometrica, 48(1), 1–48.

Author

Sameer Chawla
MSc Economics
Gokhale Institute of Politics and Economics (GIPE), Pune
Expected graduation: 2027

This repository is an academic econometrics project. Its empirical results should be interpreted as reduced-form evidence, not as a definitive structural causal model of the U.S. labour market.
