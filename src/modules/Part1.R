# PART 1: Exploratory analysis and stationarity transformation for the Building Construction IPI series
# ################################### Question 1 ################################### 
# Description of the series:
#    - "Indice CVS-CJO de la production industrielle (base 100 en 2021) - Construction de bâtiments"
#    - Monthly frequency, seasonally and working-day adjusted (CVS-CJO) volume index (NACE div. 41).
#    - Base year 2021 = 100. Values >100 indicate production above the 2021 average.

# 2. Load required packages for all the project
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("tsibble", quietly = TRUE)) install.packages("tsibble")
if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tseries", quietly = TRUE)) install.packages("tseries")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("urca", quietly = TRUE)) install.packages("urca")
if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
if (!requireNamespace("forecast", quietly = TRUE)) install.packages("forecast")
if (!requireNamespace("ellipse", quietly = TRUE)) install.packages("ellipse")


library(readr)
library(tsibble)
library(lubridate)
library(ggplot2)
library(tseries)
library(dplyr)
library(urca)
library(zoo)
library(forecast)
library(ellipse)

# 3. Import the series from a local CSV file downloaded from INSEE


data_raw <- read.csv("/Users/elouan/Desktop/Time Series project /serie_010768718_03052025/valeurs_mensuelles.csv", sep=";")
data_ts <- data_raw %>%
  # 1) rename the CSV’s ‘Date’ column to a safe name
  rename(obs_date = Date) %>%       
  
  # 2) create a true Date (we assume obs_date is "YYYY-MM")
  mutate(
    date = as.Date(paste0(obs_date, "-01"))
  ) %>%
  
  # 3) build the year‐month index
  mutate(
    ym = yearmonth(date)
  ) %>%
  
  # 4) finally turn into a tsibble
  as_tsibble(index = ym)




################################### Question 2&3 ################################### 



## 4. Analyze the raw series to check stationarity
ggplot(data_ts, aes(x = ym, y = Index)) +   
  geom_line() +                            
  labs(
    title = "Raw IPI – Building Construction",
    x     = "Date",
    y     = "Index (Base 100 = 2021)"
  ) +
  theme_minimal() 



data_ts %>%
  mutate(
    roll_mean  = zoo::rollmean(Index, 24, fill = NA),
    roll_sd    = zoo::rollapply(Index, 24, sd, fill = NA)
  ) %>%
  ggplot(aes(x = ym)) +
  geom_line(aes(y = Index)) +
  geom_line(aes(y = roll_mean),   color = "blue") +
  geom_line(aes(y = roll_sd * 10), color = "red") 

y <- data_ts$Index   # or data_ts$log_index, or data_ts$dlog_index

# Plot ACF up to lag 36 (3 years of monthly data)
acf(
  y,
  lag.max = 200,            # maximum lag to display
  main    = "ACF of IPI Series",
  ylab    = "Autocorrelation",
  xlab    = "Lag (months)"
)
# It is evident that there is a dertministic trend downward we can test it :

# add a simple time counter 1,2,…,T
data_trend <- data_ts %>%
  mutate(t = row_number())

# regress index on t
fit_lin <- lm(Index ~ t, data = data_trend)
summary(fit_lin)

#The p-value of t is very low then t is statistically significant and there is well a dertiministic trend downward



#we continue to check it with ADF test


# 1a) extract your series as a plain vector
#    (here I assume you've stored the level in data_ts$index)
y <- data_ts$Index  

#1b) We have to find the right lag order for the ADF test 
Qtests <- function(series, k, fitdf = 0) {
  pvals <- sapply(1:k, function(l) {
    if (l <= fitdf) {
      return(NA)
    } else {
      return(Box.test(series, lag = l, type = "Ljung-Box", fitdf = fitdf)$p.value)
    }
  })
  out <- data.frame(lag = 1:k, pval = pvals)
  return(out)
}

adf_until_white <- function(series, kmax = 24, type = "ct") {
  k <- 0
  noautocorr <- FALSE
  while (!noautocorr && k <= kmax) {
    cat(paste0("ADF with ", k, " lags: residuals OK? "))
    
    adf <- ur.df(series, lags = k, type = ifelse(type == "c", "drift", ifelse(type == "ct", "trend", "none")))
    
    # Run Ljung-Box tests up to lag 24, adjusting for number of regressors
    pvals_df <- Qtests(adf@testreg$residuals, 24, fitdf = length(adf@testreg$coefficients))
    
    if (sum(pvals_df$pval < 0.05, na.rm = TRUE) == 0) {
      noautocorr <- TRUE
      cat("OK ✅\n")
    } else {
      cat("nope ❌\n")
      k <- k + 1
    }
  }
  
  if (k > kmax) {
    warning("Maximum number of lags reached without achieving white residuals.")
  }
  
  return(adf)
}

adf_result <- adf_until_white(y, 24, "ct")


#After adding 3 lags, the residuals were no longer autocorrelated (based on Ljung-Box tests up to lag 24).
#We have now to perform an ADF test of laf order 3.

#ADF with 0 lags: residuals OK? nope ❌
#ADF with 1 lags: residuals OK? nope ❌
#ADF with 2 lags: residuals OK? nope ❌
#ADF with 3 lags: residuals OK? OK ✅

# 1c) run the ADF test allowing for an intercept + trend
adf_trend <- ur.df(
  y,
  type       = "trend",    # intercept + time trend
  lags = 3         # select 3 as lag order
)

# 1d) inspect the output
summary(adf_trend)


############################################### 
# Augmented Dickey-Fuller Test Unit Root Test # 
############################################### 

# Test regression trend 
# 
# 
# Call:
#   lm(formula = z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -44.821  -1.690  -0.028   2.428  20.562 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 41.17577    8.75236   4.705 4.37e-06 ***
#   z.lag.1     -0.26529    0.05629  -4.713 4.21e-06 ***
#   tt          -0.07155    0.01552  -4.611 6.63e-06 ***
#   z.diff.lag1 -0.01417    0.07098  -0.200   0.8419    
#   z.diff.lag2 -0.15866    0.06585  -2.409   0.0168 *  
#   z.diff.lag3 -0.07746    0.06458  -1.200   0.2315    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.057 on 232 degrees of freedom
# Multiple R-squared:  0.1839,	Adjusted R-squared:  0.1664 
# F-statistic: 10.46 on 5 and 232 DF,  p-value: 4.582e-09
# 
# 
# Value of test-statistic is: -4.7128 7.7748 11.2998 

# Critical values for test statistics: 
#   1pct  5pct 10pct
# tau3 -3.99 -3.43 -3.13
# phi2  6.22  4.75  4.07
# phi3  8.43  6.49  5.47



#then we can conclude we have stationnartiy but we have a drift and a trend are present !

## 1 we remove the trend by  Linear detrending

# 1a) Create a simple time index
data_ts <- data_ts %>% mutate(t = row_number())

# 1b) Fit the linear trend
fit_trend <- lm(Index ~ t, data = data_ts)

# 1c) Subtract the fitted trend
data_ts <- data_ts %>%
  mutate(
    trend_hat  = predict(fit_trend),
    detrended  = Index - trend_hat
  )

y2 <- data_ts$detrended  


#"nc" = no constant, no trend

#"c" = constant only (intercept)

#"ct" = constant + trend

#1d) We have to find the right lag order for the ADF test
adf_result <- adf_until_white(y2, 24, "c")

# ADF with 1 lags: residuals OK? nope ❌
# ADF with 2 lags: residuals OK? nope ❌
# ADF with 3 lags: residuals OK? nope ❌
# ADF with 4 lags: residuals OK? OK ✅

# 1e) run the ADF test allowing for an intercept + trend
adf_trend2 <- ur.df(
  y2,
  type = "drift",    # intercept only
  lags = 4       # choose optimal lag length
)

# 1f) inspect the output
summary(adf_trend2)


############################################### 
# Augmented Dickey-Fuller Test Unit Root Test # 
############################################### 

# Test regression drift 
# 
# 
# Call:
#   lm(formula = z.diff ~ z.lag.1 + 1 + z.diff.lag)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -45.544  -1.720   0.126   2.412  19.456 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.12848    0.39334   0.327  0.74424    
#   z.lag.1     -0.24511    0.05891  -4.161 4.48e-05 ***
#   z.diff.lag1 -0.03992    0.07410  -0.539  0.59064    
#   z.diff.lag2 -0.19219    0.07090  -2.711  0.00721 ** 
#   z.diff.lag3 -0.09978    0.06662  -1.498  0.13560    
#   z.diff.lag4 -0.08958    0.06471  -1.384  0.16763    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.051 on 231 degrees of freedom
# Multiple R-squared:  0.1886,	Adjusted R-squared:  0.171 
# F-statistic: 10.74 on 5 and 231 DF,  p-value: 2.695e-09
# 
# 
# Value of test-statistic is: -4.1607 8.6753 
# 
# Critical values for test statistics: 
#   1pct  5pct 10pct
# tau2 -3.46 -2.88 -2.57
# phi1  6.52  4.63  3.81



#recheck the graph of moving average and rolling volatility 

data_ts %>%
  mutate(
    roll_mean  = zoo::rollmean(detrended, 24, fill = NA),
    roll_sd    = zoo::rollapply(detrended, 24, sd, fill = NA)
  ) %>%
  ggplot(aes(x = ym)) +
  geom_line(aes(y = detrended)) +
  geom_line(aes(y = roll_mean),   color = "blue") +
  geom_line(aes(y = roll_sd * 10), color = "red") 

#replot the ACF


# Plot ACF up to lag 200
acf(
  y2,
  lag.max = 200,            # maximum lag to display
  main    = "ACF of IPI Series",
  ylab    = "Autocorrelation",
  xlab    = "Lag (months)"
)

#we better proof of stationnarity, let's try another method

# 2 First differencing (exactly removes any linear drift)

data_ts <- data_ts %>%
  mutate(differenced = Index - lag(Index)) %>%
  na.omit()

y3 <- data_ts$differenced  

#1b) We have to find the right lag order for the ADF test

adf_result <- adf_until_white(y3, 24, "nc")

# ADF with 0 lags: residuals OK? nope ❌
# ADF with 1 lags: residuals OK? nope ❌
# ADF with 2 lags: residuals OK? nope ❌
# ADF with 3 lags: residuals OK? nope ❌
# ADF with 4 lags: residuals OK? nope ❌
# ADF with 5 lags: residuals OK? OK ✅

# 1c) run the ADF test allowing for an intercept + trend
adf_trend3 <- ur.df(
  y3,
  type       = "none",    # no intercept and no trend
  lags = 5       # choose optimal lag length
)

# 1d) inspect the output
summary(adf_trend3)


############################################### 
# Augmented Dickey-Fuller Test Unit Root Test # 
############################################### 

# Test regression none 
# 
# 
# Call:
#   lm(formula = z.diff ~ z.lag.1 - 1 + z.diff.lag)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -46.075  -2.116  -0.391   1.803  21.905 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# z.lag.1     -2.28954    0.24724  -9.260  < 2e-16 ***
#   z.diff.lag1  1.05822    0.21746   4.866 2.12e-06 ***
#   z.diff.lag2  0.69133    0.18348   3.768 0.000209 ***
#   z.diff.lag3  0.43989    0.14452   3.044 0.002609 ** 
#   z.diff.lag4  0.22703    0.10324   2.199 0.028883 *  
#   z.diff.lag5  0.09011    0.06558   1.374 0.170791    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.24 on 229 degrees of freedom
# Multiple R-squared:  0.6135,	Adjusted R-squared:  0.6034 
# F-statistic: 60.59 on 6 and 229 DF,  p-value: < 2.2e-16
# 
# 
# Value of test-statistic is: -9.2604 
# 
# Critical values for test statistics: 
#   1pct  5pct 10pct
# tau1 -2.58 -1.95 -1.62


#This time, we observe no more trend and drift, and the serie is stationary without a trend and a constant ! 


#recheck the graph of moving average and rolling volatility 

data_ts %>%
  mutate(
    roll_mean  = zoo::rollmean(differenced, 24, fill = NA),
    roll_sd    = zoo::rollapply(differenced, 24, sd, fill = NA)
  ) %>%
  ggplot(aes(x = ym)) +
  geom_line(aes(y = differenced)) +
  geom_line(aes(y = roll_mean),   color = "blue") +
  geom_line(aes(y = roll_sd * 10), color = "red") 


# We have a better stationarity around the mean, but it is a bit worse for the volatility


#Replot the ACF


# Plot ACF up to lag 200
acf(
  y3,
  lag.max = 200,            # maximum lag to display
  main    = "ACF of IPI Series",
  ylab    = "Autocorrelation",
  xlab    = "Lag (months)"
)

#We have clearly better results for the ACF function


ggplot(data_ts, aes(x = ym, y = differenced)) +   
  geom_line() +                            
  labs(
    title = "Transformed IPI – Building Construction",
    x     = "Date",
    y     = "Index (Base 100 = 2021)"
  ) +
  theme_minimal() 

# We observe that the new series does not have a drift and a trend anymore.


# End of Part 1
