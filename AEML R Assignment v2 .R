# =============================================================================
# HIGH FREQUENCY REALIZED VOLATILITY MODELS ASSINMENT
# Bukho Ndziweni, Samson Mashabane

# Data set used: Boeing (BA) Intraday data: Jan 2003 to Dec 2008
# Observations: 605671, at 1-minute frequency
# Columns in data set: TICKER, DATE, TIME, OPEN, HIGH, LOW, CLOSE, VOLUME
# =============================================================================


# PACKAGE INSTALLATION 
install.packages(c("highfrequency", "xts", "zoo", "rugarch",
                   "ggplot2", "moments", "nortest", "data.table"))

# Loading all packages
library(highfrequency)  # Realized volatility estimators
library(xts)            # Extensible time series objects
library(zoo)            # Rolling window utilities
library(rugarch)        # GARCH model fitting
library(ggplot2)        # Plotting
library(moments)        # Jarque-Bera normality test
library(nortest)        # Lilliefors normality test
library(data.table)     # Fast CSV reading (fread)


# =============================================================================
# QUESTION 1: Loading Boeing intraday stock price data
# =============================================================================

# Since fread() is much faster than read.csv() for large files (~600k rows), we use that
ba_raw <- fread(
  "BA.asc",            
  header = TRUE,       
  sep    = ","         
)

# We preview to confirm correct loading
cat("Dimensions:", nrow(ba_raw), "rows x", ncol(ba_raw), "columns\n")
cat("Column names:", paste(names(ba_raw), collapse = ", "), "\n")
print(head(ba_raw, 5))


# =============================================================================
# QUESTION 2: Creation of an xts object for the highfrequency package
# =============================================================================

# Here, we combine DATE and TIME into a single POSIXct timestamp
# The DATE format is MM/DD/YYYY, TIME is HH:MM:SS
ba_datetime <- as.POSIXct(
  paste(ba_raw$DATE, ba_raw$TIME),  
  format = "%m/%d/%Y %H:%M:%S",    
  tz     = "America/New_York" 
)

# Next, we build the price xts using CLOSE prices -- used for return-based estimators
ba_price <- xts(
  ba_raw$CLOSE,      
  order.by = ba_datetime 
)

# Here we build OHLCV xts because it is needed for the Realized Range estimator (Q5d)
ba_ohlcv <- xts(
  ba_raw[, .(OPEN, HIGH, LOW, CLOSE, VOLUME)], 
  order.by = ba_datetime
)

cat("xts price object created:", nrow(ba_price), "observations\n")
cat("Date range:", format(start(ba_price)), "to", format(end(ba_price)), "\n")


# =============================================================================
# QUESTION 3: Dropping weekend days
# =============================================================================

# .indexwday: 0 = Sunday, 1 = Monday, ..., 5 = Friday, 6 = Saturday
# Keeping only Monday through Friday (weekdays 1 to 5)
ba_price  <- ba_price[.indexwday(ba_price)  %in% 1:5]
ba_ohlcv  <- ba_ohlcv[.indexwday(ba_ohlcv) %in% 1:5]

cat("After dropping weekends:", nrow(ba_price), "observations remain\n")


# =============================================================================
# QUESTION 4: The visual inspection and removal of anomalous observations
# =============================================================================

# ---- 4a: The Plot of the raw intraday price series ----------------------------------
plot(ba_price,
     main = "Boeing (BA) -- Raw 1-Minute Close Prices (2003-2008)",
     ylab = "Price (USD)",
     col  = "steelblue",
     lwd  = 0.3)

# ---- 4b: Computation of 1-minute log returns ----------------------------------------
# Log return at time t: r_t = ln(P_t) - ln(P_{t-1})
ba_ret_raw <- diff(log(ba_price))  
ba_ret_raw <- na.omit(ba_ret_raw) 

# ---- 4c: Detection and removal of abnormal returns ----------------------------------
# Following Brownlees & Gallo (2008): flag returns whose absolute value exceeds
# a multiple of the median absolute deviation (MAD) -- a robust threshold.
# We use 25 * MAD, a common choice in the microstructure literature.

mad_scale  <- mad(as.vector(ba_ret_raw), na.rm = TRUE)  
threshold  <- 25 * mad_scale  

cat("\nReturn threshold (25 * MAD):", round(threshold, 6), "\n")

# Identifying observations exceeding the threshold
outlier_flag <- abs(as.vector(ba_ret_raw)) > threshold
n_outliers   <- sum(outlier_flag, na.rm = TRUE)
cat("Abnormal return observations found:", n_outliers, "\n")

if (n_outliers > 0) {
  # Getting time stamps of outlier returns (these correspond to the *second* price
  # in each return pair, i.e. the price that caused the spike)
  bad_times <- index(ba_ret_raw)[outlier_flag]
  print(data.frame(
    Timestamp = bad_times,
    Return    = as.vector(ba_ret_raw)[outlier_flag]
  ))
  # Removing those timestamps from the price and OHLCV series
  ba_price_clean <- ba_price[!(index(ba_price) %in% bad_times)]
  ba_ohlcv_clean <- ba_ohlcv[!(index(ba_ohlcv) %in% bad_times)]
} else {
  cat("No anomalies found -- using full dataset\n")
  ba_price_clean <- ba_price
  ba_ohlcv_clean <- ba_ohlcv
}

# ---- 4d: Re-ploting a cleaned price series ----------------------------------------
plot(ba_price_clean,
     main = "Boeing (BA) -- Cleaned 1-Minute Close Prices (2003-2008)",
     ylab = "Price (USD)",
     col  = "darkgreen",
     lwd  = 0.3)


# =============================================================================
# QUESTION 5: Realized Volatility Estimators
# =============================================================================

# ---- 5a: Classical RV -- Andersen & Bollerslev (1998) -----------------------
rv_classical <- rCov(
  rData       = ba_price_clean,
  alignBy     = "minutes",
  alignPeriod = 5,
  makeReturns = TRUE
)
colnames(rv_classical) <- "RV_Classical"
cat("Classical RV -- first 5 trading days:\n")
print(head(rv_classical, 5))

# ---- 5b: Two-Time-Scale Estimator (loop over days, pass prices) -------------
trading_dates <- unique(as.Date(index(ba_price_clean)))
cat("Number of trading days:", length(trading_dates), "\n")

tts_values <- numeric(length(trading_dates))

for (i in seq_along(trading_dates)) {
  day_price <- ba_price_clean[as.Date(index(ba_price_clean)) == trading_dates[i]]
  if (nrow(day_price) > 10) {
    tryCatch({
      tts_values[i] <- as.numeric(rTSCov(pData = day_price, K = 20, J = 1))
    }, error = function(e) {
      tts_values[i] <<- NA
    })
  } else {
    tts_values[i] <- NA
  }
  if (i %% 100 == 0) cat("Processed", i, "of", length(trading_dates), "days\n")
}

rv_tts <- xts(tts_values, order.by = trading_dates)
colnames(rv_tts) <- "RV_TTS"
cat("Two-Time-Scale RV -- first 5 trading days:\n")
print(head(rv_tts, 5))

# ---- 5c: Realized Kernel Estimator ------------------------------------------
rv_kernel <- rKernelCov(
  rData       = ba_price_clean,
  alignBy     = "minutes",
  alignPeriod = 1,
  makeReturns = TRUE,
  kernelType  = "parzen"
)
colnames(rv_kernel) <- "RV_Kernel"
cat("Realized Kernel RV -- first 5 trading days:\n")
print(head(rv_kernel, 5))

# ---- 5d: Realized Range Estimator (fully manual) ----------------------------
# Here we aggregate High and Low to 5-minute bars using endpoints() and period.apply()

# Extract High and Low columns
ba_high <- ba_ohlcv_clean[, "HIGH"]
ba_low  <- ba_ohlcv_clean[, "LOW"]

# 5-minute endpoints (breakpoints every 5 minutes)
ep <- endpoints(ba_high, on = "minutes", k = 5)

# Applying max to High and min to Low within each 5-minute window
ba_high_5min <- period.apply(ba_high, INDEX = ep, FUN = max)
ba_low_5min  <- period.apply(ba_low,  INDEX = ep, FUN = min)

cat("High 5-min rows:", nrow(ba_high_5min), "\n")
cat("Low 5-min rows:",  nrow(ba_low_5min),  "\n")
print(head(ba_high_5min, 3))

# Parkinson (1980) scaled squared log range for each 5-minute bar
log_range_sq        <- (log(ba_high_5min) - log(ba_low_5min))^2
log_range_sq_scaled <- log_range_sq / (4 * log(2))

# Sum within each calendar day
rv_range <- apply.daily(log_range_sq_scaled, sum, na.rm = TRUE)
colnames(rv_range) <- "RV_Range"
cat("Realized Range RV -- first 5 trading days:\n")
print(head(rv_range, 5))

# ---- Ploting all four estimators ------------------------------------------------
common_rv_dates <- Reduce(
  intersect,
  list(as.character(as.Date(index(rv_classical))),
       as.character(as.Date(index(rv_tts))),
       as.character(as.Date(index(rv_kernel))),
       as.character(as.Date(index(rv_range))))
)

common_rv_dates <- as.Date(common_rv_dates)
cat("Common trading days:", length(common_rv_dates), "\n")

rv_classical_v <- as.vector(rv_classical[as.Date(index(rv_classical)) %in% common_rv_dates])
rv_tts_v       <- as.vector(rv_tts[as.Date(index(rv_tts)) %in% common_rv_dates])
rv_kernel_v    <- as.vector(rv_kernel[as.Date(index(rv_kernel)) %in% common_rv_dates])
rv_range_v     <- as.vector(rv_range[as.Date(index(rv_range)) %in% common_rv_dates])

rv_df_long <- data.frame(
  Date      = rep(common_rv_dates, 4),
  RV        = c(rv_classical_v, rv_tts_v, rv_kernel_v, rv_range_v),
  Estimator = rep(c("Classical", "TTS", "Kernel", "Range"),
                  each = length(common_rv_dates))
)

ggplot(rv_df_long, aes(x = Date, y = RV, color = Estimator)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  facet_wrap(~Estimator, scales = "free_y", ncol = 2) +
  labs(title = "Boeing Daily Realized Variance: Four Estimators (2003-2008)",
       x = "Date", y = "Realized Variance") +
  theme_minimal() +
  theme(legend.position = "none")

# =============================================================================
# QUESTION 6: Normality of daily returns standardized by RV
# =============================================================================

# ---- Q6: Daily returns standardized by RV -- normality test -----------------
daily_prices <- aggregateTS(
  ba_price_clean,
  alignBy     = "minutes",
  alignPeriod = 390
)

daily_ret <- diff(log(daily_prices))
daily_ret <- na.omit(daily_ret)

common_dates_q6 <- intersect(
  as.character(as.Date(index(daily_ret))),
  as.character(as.Date(index(rv_classical)))
)
common_dates_q6 <- as.Date(common_dates_q6)

ret_q6 <- as.vector(daily_ret[as.Date(index(daily_ret)) %in% common_dates_q6])
rv_q6  <- as.vector(rv_classical[as.Date(index(rv_classical)) %in% common_dates_q6])

valid  <- rv_q6 > 0 & !is.na(rv_q6)
ret_valid <- ret_q6[valid]
rv_valid  <- rv_q6[valid]
cat("Valid trading days:", sum(valid), "\n")

std_ret <- ret_valid / sqrt(rv_valid)
cat("Skewness:", round(skewness(std_ret), 4), "\n")
cat("Excess kurtosis:", round(kurtosis(std_ret) - 3, 4), "\n")

jb_result <- jarque.test(std_ret)
print(jb_result)

ll_result <- lillie.test(std_ret)
print(ll_result)

qqnorm(std_ret, main = "QQ-Plot: Daily Returns Standardized by sqrt(RV)",
       pch = 16, cex = 0.5, col = "steelblue")
qqline(std_ret, col = "red", lwd = 2)


# =============================================================================
# QUESTION 7: Jump test and visualization
# =============================================================================

# ---- 7a: BNS Jump Test (Barndorff-Nielsen & Shephard, 2004/2006) ------------
# The strategy: compare Realized Variance (RV) to Bipower Variation (BV).
# BV is robust to jumps (uses products of adjacent absolute returns).
# A significantly large RV relative to BV signals a jump.
#
# Test statistic: z_t = (RV_t - BV_t) / (sqrt(variance_ratio * BV_t^2 / n))
# We use the ratio statistic: J_t = 1 - BV_t / RV_t
# Under no jump: J_t ~ 0; under a jump: J_t > 0.

# ---- Q7: Jump test and visualization ----------------------------------------
rv_bv <- rBPCov(
  rData       = ba_price_clean,
  alignBy     = "minutes",
  alignPeriod = 5,
  makeReturns = TRUE
)
colnames(rv_bv) <- "BV"

common_jump_dates <- intersect(
  as.character(as.Date(index(rv_classical))),
  as.character(as.Date(index(rv_bv)))
)
common_jump_dates <- as.Date(common_jump_dates)

rv_j <- as.vector(rv_classical[as.Date(index(rv_classical)) %in% common_jump_dates])
bv_j <- as.vector(rv_bv[as.Date(index(rv_bv)) %in% common_jump_dates])

# BNS ratio statistic
mu1       <- sqrt(2 / pi)
n_obs     <- 78
theta     <- (pi^2 / 4 + pi - 5)
rel_jump  <- pmax(1 - bv_j / rv_j, 0)
z_bns     <- sqrt(n_obs) * rel_jump / sqrt(theta * pmax(rv_j^2, 1e-12) / bv_j^2)

jump_flag <- z_bns > qnorm(0.99)
n_jumps   <- sum(jump_flag, na.rm = TRUE)
cat("Number of jump days detected:", n_jumps, "\n")

rv_jump_df <- data.frame(
  Date = common_jump_dates,
  RV   = rv_j,
  Jump = jump_flag
)

ggplot(rv_jump_df, aes(x = Date, y = RV)) +
  geom_line(color = "steelblue", linewidth = 0.5) +
  geom_point(
    data   = rv_jump_df[rv_jump_df$Jump & !is.na(rv_jump_df$Jump), ],
    aes(x = Date, y = RV),
    shape  = 1, color = "red", size = 3, stroke = 1.5
  ) +
  labs(
    title    = "Boeing Daily Realized Variance (2003-2008)",
    subtitle = "Red circles indicate jump days (BNS test, 1% level)",
    x = "Date", y = "Realized Variance"
  ) +
  theme_minimal()


# =============================================================================
# QUESTION 8: Signature plot -- RV at different sampling frequencies
# =============================================================================

# ---- Q8: Signature plot -----------------------------------------------------
freqs <- c(1, 5, 10, 20, 30, 40, 50, 60)
rv_sd <- numeric(length(freqs))

for (i in seq_along(freqs)) {
  rv_f <- rCov(
    rData       = ba_price_clean,
    alignBy     = "minutes",
    alignPeriod = freqs[i],
    makeReturns = TRUE
  )
  rv_sd[i] <- sd(as.vector(rv_f), na.rm = TRUE)
  cat("Frequency:", freqs[i], "min -- SD of RV:", round(rv_sd[i], 8), "\n")
}

sig_df <- data.frame(
  Frequency_min = freqs,
  SD_RV         = rv_sd
)

print(sig_df)

ggplot(sig_df, aes(x = Frequency_min, y = SD_RV)) +
  geom_line(color = "darkorange", linewidth = 1.2) +
  geom_point(color = "darkorange", size = 3) +
  scale_x_continuous(breaks = freqs) +
  labs(
    title    = "Volatility Signature Plot",
    subtitle = "SD of Daily RV inflates at high frequency due to microstructure noise",
    x        = "Sampling Frequency (minutes)",
    y        = "Sample Standard Deviation of Daily RV"
  ) +
  theme_minimal()

# Interpretation
# at 1-minute, microstructure noise (bid-ask bounce, rounding)
# inflates RV and its variability. As frequency coarsens (longer intervals),
# noise averages out and RV stabilizes -- the classic signature plot pattern.
# The "elbow" near 5 minutes is why 5-min is the industry convention.


# =============================================================================
# QUESTION 9: HAR Model calibration and in-sample forecast
# =============================================================================
# Corsi (2009) HAR model:
#   RV_t = c + beta_d * RV_{t-1}
#             + beta_w * RV^(5)_{t-1}   [5-day average]
#             + beta_m * RV^(22)_{t-1}  [22-day average]
#             + epsilon_t
#
# Captures long-memory-like persistence via three cascade components.

# We use the classical 5-minute RV series
# ---- Q9: HAR Model (fast manual implementation using lm) --------------------
rv_vec <- as.numeric(rv_classical)
n      <- length(rv_vec)

# Daily, weekly and monthly components
rv_d <- rv_vec                                  
rv_w <- filter(rv_vec, rep(1/5,  5),  sides = 1)
rv_m <- filter(rv_vec, rep(1/22, 22), sides = 1)

# Aligning dependent variable starts at t=23, regressors at t=22
y    <- rv_vec[23:n]          
x_d  <- rv_d[22:(n-1)]      
x_w  <- rv_w[22:(n-1)]       
x_m  <- rv_m[22:(n-1)]       

har_df_lm <- data.frame(y = y, x_d = x_d, x_w = x_w, x_m = x_m)
har_df_lm <- na.omit(har_df_lm)

# Fitting HAR via OLS
har_lm <- lm(y ~ x_d + x_w + x_m, data = har_df_lm)

cat("\n--- HAR Model Summary ---\n")
print(summary(har_lm))

coefs <- coef(har_lm)
cat("\nIntercept:", round(coefs[1], 8), "\n")
cat("Beta daily:", round(coefs[2], 4), "\n")
cat("Beta weekly:", round(coefs[3], 4), "\n")
cat("Beta monthly:", round(coefs[4], 4), "\n")
cat("Sum of betas (persistence):", round(sum(coefs[-1]), 4), "\n")

cat("\nBox-Ljung test on residuals:\n")
print(Box.test(residuals(har_lm), lag = 10, type = "Ljung-Box"))

# In-sample fitted values plot
har_dates_lm <- tail(as.Date(index(rv_classical)), nrow(har_df_lm))

har_plot_df <- data.frame(
  Date   = har_dates_lm,
  Actual = har_df_lm$y,
  Fitted = fitted(har_lm)
)

ggplot(har_plot_df, aes(x = Date)) +
  geom_line(aes(y = Actual, color = "Actual RV"), linewidth = 0.4) +
  geom_line(aes(y = Fitted, color = "HAR Forecast"),
            linewidth = 0.4, linetype = "dashed") +
  scale_color_manual(values = c("Actual RV"    = "steelblue",
                                "HAR Forecast" = "red")) +
  labs(title = "HAR Model: In-Sample 1-Day-Ahead RV Forecast vs Actual",
       x = "Date", y = "Realized Variance", color = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# =============================================================================
# QUESTION 10: GARCH model on daily returns
# =============================================================================
# Model choice: GJR-GARCH(1,1) with Student-t errors.
#
# Justification would be the following;
#   - GARCH(1,1): parsimonious; captures volatility clustering (well established
#     in the equity literature -- Bollerslev, 1986).
#   - GJR asymmetry term (gamma): negative shocks increase variance more than
#     positive shocks of the same magnitude (leverage effect, Nelson, 1991).
#     For equity returns this is almost always significant.
#   - Student-t distribution: daily returns exhibit heavier tails than Gaussian;
#     t-distribution with estimated degrees of freedom fits this better.

# ---- Q10: GJR-GARCH(1,1) with Student-t errors ------------------------------
daily_prices <- aggregateTS(
  ba_price_clean,
  alignBy     = "minutes",
  alignPeriod = 390
)
daily_ret    <- diff(log(daily_prices))
daily_ret    <- na.omit(daily_ret)
daily_ret_vec <- as.numeric(daily_ret)

garch_spec <- ugarchspec(
  variance.model = list(
    model      = "gjrGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder    = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

garch_fit <- ugarchfit(
  spec   = garch_spec,
  data   = daily_ret_vec,
  solver = "hybrid"
)

cat("\n--- GJR-GARCH(1,1) Summary ---\n")
print(garch_fit)

garch_std_resid <- residuals(garch_fit, standardize = TRUE)
qqnorm(garch_std_resid,
       main = "QQ-Plot: GJR-GARCH Standardized Residuals",
       pch = 16, cex = 0.5, col = "steelblue")
qqline(garch_std_resid, col = "red", lwd = 2)

# =============================================================================
# QUESTION 11: Model comparison via Mincer-Zarnowitz regression
# =============================================================================
# Mincer & Zarnowitz (1969) regression:
#   RV_t = alpha + beta * sigma^2_hat_t + u_t
#
# A perfect forecast gives alpha = 0 and beta = 1 (unbiased).
# We compare: (a) HAR forecast of RV, (b) GARCH conditional variance.
# Both are compared against the 5-minute classical RV as the "true" volatility.

# =============================================================================
# QUESTION 11: Mincer-Zarnowitz Regression
# =============================================================================

# Fixing daily returns using endpoints-based aggregation
daily_prices_fix  <- ba_price_clean[endpoints(ba_price_clean, on = "days")]
daily_ret_fix     <- diff(log(daily_prices_fix))
daily_ret_fix     <- na.omit(daily_ret_fix)
daily_ret_vec_fix <- as.numeric(daily_ret_fix)
garch_dates_fix   <- as.Date(index(daily_ret_fix))
cat("Daily returns:", length(daily_ret_vec_fix), "\n")

# Refit GARCH on correct daily returns
garch_spec <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
garch_fit <- ugarchfit(
  spec   = garch_spec,
  data   = daily_ret_vec_fix,
  solver = "hybrid"
)
cat("\n--- GJR-GARCH Summary ---\n")
print(garch_fit)

# GARCH QQ plot
garch_std_resid <- residuals(garch_fit, standardize = TRUE)
qqnorm(garch_std_resid,
       main = "QQ-Plot: GJR-GARCH Standardized Residuals",
       pch = 16, cex = 0.5, col = "steelblue")
qqline(garch_std_resid, col = "red", lwd = 2)

# Extract GARCH conditional variance
garch_var_fix <- as.numeric(sigma(garch_fit))^2

# HAR forecasts and dates
har_dates_mz  <- tail(as.Date(index(rv_classical)), nrow(har_df_lm))
har_fcst_mz   <- as.numeric(fitted(har_lm))
rv_true_dates <- as.Date(index(rv_classical))
rv_true       <- as.numeric(rv_classical)

# Find common dates across all three
mz_common <- Reduce(
  intersect,
  list(as.character(har_dates_mz),
       as.character(garch_dates_fix),
       as.character(rv_true_dates))
)
mz_common <- as.Date(mz_common)
cat("Common MZ dates:", length(mz_common), "\n")

rv_mz  <- rv_true[rv_true_dates %in% mz_common]
har_mz <- har_fcst_mz[har_dates_mz %in% mz_common]
gch_mz <- garch_var_fix[garch_dates_fix %in% mz_common]

cat("Lengths -- rv:", length(rv_mz),
    "har:", length(har_mz),
    "garch:", length(gch_mz), "\n")

# MZ regressions
mz_har   <- lm(rv_mz ~ har_mz)
mz_garch <- lm(rv_mz ~ gch_mz)

cat("\n--- MZ: HAR vs 5-min RV ---\n")
print(summary(mz_har))

cat("\n--- MZ: GARCH vs 5-min RV ---\n")
print(summary(mz_garch))

cat("\n--- Comparison Table ---\n")
mz_compare <- data.frame(
  Model = c("HAR", "GJR-GARCH"),
  Alpha = round(c(coef(mz_har)[1],   coef(mz_garch)[1]),   6),
  Beta  = round(c(coef(mz_har)[2],   coef(mz_garch)[2]),   4),
  R2    = round(c(summary(mz_har)$r.squared,
                  summary(mz_garch)$r.squared), 4)
)
print(mz_compare)

# Scatter plot
mz_plot_df <- rbind(
  data.frame(Model = "HAR",       Forecast = har_mz, Actual = rv_mz),
  data.frame(Model = "GJR-GARCH", Forecast = gch_mz, Actual = rv_mz)
)

ggplot(mz_plot_df, aes(x = Forecast, y = Actual)) +
  geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "red") +
  facet_wrap(~Model, scales = "free_x") +
  labs(
    title    = "Mincer-Zarnowitz: Forecast vs Actual Realized Variance",
    subtitle = "Dashed line = perfect forecast (alpha=0, beta=1)",
    x = "Forecast", y = "Actual RV (5-min)"
  ) +
  theme_minimal()
# =============================================================================
# QUESTION 12, The bonus question: Bitcoin and Ethereum realized correlation
# =============================================================================
# Uses the binancer package to pull 30-min Binance data from Jan 2025 to May 2026.

# install.packages("binancer")
library(binancer) 

# ---- 12a: We sart byfetching 30-minute BTC and ETH price data ---------------
# binance_klines() returns OHLCV kline (candlestick) data from Binance
# Binance allows up to 1000 rows per call; we loop over the date range.
# Date range: Jan 1, 2025 to May 1, 2026

# ---- Q12: Bitcoin and Ethereum realized correlation -------------------------
# Helper function to fetch data in chunks (Binance limits per call)
fetch_binance <- function(symbol, start, end, interval = "30m") {
  all_data <- list()
  current  <- as.POSIXct(start, tz = "UTC")
  end_time <- as.POSIXct(end,   tz = "UTC")
  
  while (current < end_time) {
    chunk <- binance_klines(
      symbol     = symbol,
      interval   = interval,
      start_time = current,
      end_time   = end_time,
      limit      = 1000
    )
    if (is.null(chunk) || nrow(chunk) == 0) break
    all_data  <- c(all_data, list(chunk))
    current   <- max(chunk$open_time) + 1
    cat("Fetched up to:", format(max(chunk$open_time)), "\n")
  }
  rbindlist(all_data)
}

cat("Fetching BTC data...\n")
btc_raw <- fetch_binance("BTCUSDT", "2025-01-01", "2026-05-01")
cat("BTC rows:", nrow(btc_raw), "\n")

cat("Fetching ETH data...\n")
eth_raw <- fetch_binance("ETHUSDT", "2025-01-01", "2026-05-01")
cat("ETH rows:", nrow(eth_raw), "\n")

# Building xts objects from BTC and ETH close prices
btc_price_xts <- xts(as.numeric(btc_raw$close),
                     order.by = as.POSIXct(btc_raw$open_time, tz = "UTC"))
eth_price_xts <- xts(as.numeric(eth_raw$close),
                     order.by = as.POSIXct(eth_raw$open_time, tz = "UTC"))

# Computing 30-minute log returns
btc_ret <- diff(log(btc_price_xts)); btc_ret <- na.omit(btc_ret)
eth_ret <- diff(log(eth_price_xts)); eth_ret <- na.omit(eth_ret)

# Debuging: check timezone and class
cat("BTC index class:", class(index(btc_ret)), "\n")
cat("ETH index class:", class(index(eth_ret)), "\n")
cat("BTC index timezone:", attr(index(btc_ret), "tzone"), "\n")
cat("ETH index timezone:", attr(index(eth_ret), "tzone"), "\n")

# Align to common timestamps using merge (safer for xts)
aligned_data <- merge(btc_ret, eth_ret, join = "inner")
btc_aligned <- aligned_data[, 1]
eth_aligned <- aligned_data[, 2]
common_crypto <- index(btc_aligned)

cat("Common 30-min observations:", length(common_crypto), "\n")

# Daily realized variances and covariance
rv_btc  <- apply.daily(btc_aligned^2,              sum, na.rm = TRUE)
rv_eth  <- apply.daily(eth_aligned^2,              sum, na.rm = TRUE)
rcov    <- apply.daily(btc_aligned * eth_aligned,  sum, na.rm = TRUE)

# Align to common daily dates using merge (safer for xts)
aligned_daily <- merge(rv_btc, rv_eth, rcov, join = "inner")
rv_btc_v <- as.vector(aligned_daily[, 1])
rv_eth_v <- as.vector(aligned_daily[, 2])
rcov_v   <- as.vector(aligned_daily[, 3])
common_daily <- index(aligned_daily)

# Realized correlation
rcor_v <- rcov_v / sqrt(rv_btc_v * rv_eth_v)
rcor_df <- data.frame(
  Date = as.Date(common_daily),
  RCor = rcor_v
)

cat("\nSummary of BTC-ETH realized correlation:\n")
print(summary(rcor_v))

# Plot
ggplot(rcor_df, aes(x = Date, y = RCor)) +
  geom_line(color = "purple", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", span = 0.3,
              color = "orange", se = FALSE, linewidth = 1) +
  labs(
    title    = "Realized Correlation: ETH vs BTC (30-min, Jan 2025 to May 2026)",
    subtitle = "Orange line = LOESS trend",
    x = "Date", y = "Realized Correlation"
  ) +
  ylim(-1, 1) +
  theme_minimal()



