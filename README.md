# High-Frequency-Realized-Volstility-Models

# High Frequency Realized Volatility Models
### Applied Econometrics and Machine Learning | MFE 2026 | University of Johannesburg

---

## Overview

This repository contains the full empirical analysis for the High Frequency Realized Volatility Models assignment, submitted in partial fulfilment of the Master of Financial Engineering (MFE) at the University of Johannesburg.

The analysis applies realized volatility estimation techniques to one-minute intraday price data for the Boeing Company (BA), traded on the NYSE from January 2003 to December 2008 (605,670 observations). A supplementary analysis examines the realized correlation between Bitcoin and Ethereum using 30-minute data from January 2025 to May 2026.

---

## Repository Structure

---

## Data

| Field | Description |
|---|---|
| Source | Provided by course instructor |
| Asset | Boeing Company (BA), NYSE |
| Frequency | 1-minute OHLCV |
| Period | January 2, 2003 to December 31, 2008 |
| Observations | 605,670 (raw), 1,553 trading days |
| Columns | TICKER, DATE, TIME, OPEN, HIGH, LOW, CLOSE, VOLUME |

---

## Analysis Summary

| Question | Topic |
|---|---|
| Q1-Q2 | Data loading and xts object construction |
| Q3 | Weekend removal |
| Q4 | Data cleaning via 25-times-MAD filter (Brownlees and Gallo, 2006) |
| Q5 | Four RV estimators: Classical, Two-Time-Scale, Realized Kernel, Realized Range |
| Q6 | Normality test on RV-standardised daily returns (Jarque-Bera, Lilliefors) |
| Q7 | Jump detection via BNS test (Barndorff-Nielsen and Shephard, 2006) |
| Q8 | Volatility signature plot across 1 to 60 minute sampling frequencies |
| Q9 | HAR model calibration and in-sample forecast (Corsi, 2009) |
| Q10 | GJR-GARCH(1,1) with Student-t errors (Glosten, Jagannathan and Runkle, 1993) |
| Q11 | Mincer-Zarnowitz model comparison: HAR vs GJR-GARCH |
| Q12 | BTC-ETH realized correlation using Binance 30-minute data (2025-2026) |

---

## Key Results

- All four RV estimators show consistent dynamics with a pronounced spike during the 2008 Global Financial Crisis
- Standardised daily returns reject normality (Jarque-Bera and Lilliefors, p < 0.01), consistent with the presence of price jumps
- The BNS jump test identifies numerous jump days concentrated in 2003 to 2004 and late 2008
- HAR outperforms GJR-GARCH in the Mincer-Zarnowitz comparison (R-squared: 0.636 vs 0.576)
- BTC-ETH realized correlation averages above 0.80 over 2025 to 2026 with a gradual upward trend

---

## Requirements

**R version:** 4.6.0 or later

**R packages:**
```r
install.packages(c(
  "highfrequency",
  "xts",
  "zoo",
  "rugarch",
  "ggplot2",
  "moments",
  "nortest",
  "data.table",
  "binancer"   # bonus question only
))
```

---

## How to Run

1. Clone the repository
2. Place `BA.asc` in the `data/` folder
3. Open `code/hf_rv_assignment.R` in RStudio
4. Set your working directory to the `data/` folder:
```r
setwd("path/to/HF-Realized-Volatility/data")
```
5. Run the script section by section from Q1 to Q12

---

## References

- Andersen, T.G. and Bollerslev, T. (1998) International Economic Review, 39(4), pp. 885-905
- Barndorff-Nielsen, O.E. and Shephard, N. (2006) Journal of Financial Econometrics, 4(1), pp. 1-30
- Corsi, F. (2009) Journal of Financial Econometrics, 7(2), pp. 174-196
- Glosten, L.R., Jagannathan, R. and Runkle, D.E. (1993) Journal of Finance, 48(5), pp. 1779-1801
- Zhang, L., Mykland, P.A. and Ait-Sahalia, Y. (2005) Journal of the American Statistical Association, 100(472), pp. 1394-1411

---

## Author

[Bukho Ndziweni and Samson Mashabane] | MFE Student | University of Johannesburg | 2026
