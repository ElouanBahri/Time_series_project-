# PART 1: Exploratory analysis and stationarity transformation for the Building Construction IPI series
# ################################### Question 1 ################################### 
# Description of the series:
#    - "Indice CVS-CJO de la production industrielle (base 100 en 2021) - Construction de bâtiments"
#    - Monthly frequency, seasonally and working-day adjusted (CVS-CJO) volume index (NACE div. 41).
#    - Base year 2021 = 100. Values >100 indicate production above the 2021 average.
#    - Because the index is multiplicative in nature, log-transformation could help to stabilize variance.

# 2. Load required packages
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("tsibble", quietly = TRUE)) install.packages("tsibble")
if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tseries", quietly = TRUE)) install.packages("tseries")

library(readr)
library(tsibble)
library(lubridate)
library(ggplot2)
library(tseries)
library(dplyr)
install.packages("urca")
library(urca)
library(zoo)


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

# Test for unit root on the series
adf_log <- adf.test(data_ts$Index, alternative = "stationary")
print(adf_log)  # If p-value > 0.05, The series is non-stationary


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

# 1b) run the ADF test allowing for an intercept + trend
adf_trend <- ur.df(
  y,
  type       = "trend",    # intercept + time trend
  selectlags = "AIC"       # choose optimal lag length
)

# 1c) inspect the output
summary(adf_trend)


############################################### 
# Augmented Dickey-Fuller Test Unit Root Test # 
############################################### 

#Test regression trend 


# Call:
#   lm(formula = z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -44.050  -1.712  -0.233   2.581  26.117 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 49.45258    7.63379   6.478 5.35e-10 ***
#   z.lag.1     -0.31959    0.04922  -6.493 4.91e-10 ***
#   tt          -0.08428    0.01375  -6.129 3.67e-09 ***
#   z.diff.lag   0.05077    0.06417   0.791     0.43    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.093 on 236 degrees of freedom
# Multiple R-squared:  0.1621,	Adjusted R-squared:  0.1514 
# F-statistic: 15.21 on 3 and 236 DF,  p-value: 4.394e-09
# 
# 
# Value of test-statistic is: -6.4934 14.239 21.2186 
# 
# Critical values for test statistics: 
#   1pct  5pct 10pct
# tau3 -3.99 -3.43 -3.13
# phi2  6.22  4.75  4.07
# phi3  8.43  6.49  5.47


#then we can conclude that a AR(1) could fit and that we have stationnartiy but we have a drift and a trend are present !

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

# 1b) run the ADF test allowing for an intercept + trend
adf_trend2 <- ur.df(
  y2,
  type       = "trend",    # intercept + time trend
  selectlags = "AIC"       # choose optimal lag length
)

# 1c) inspect the output
summary(adf_trend2)


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
# -44.050  -1.712  -0.233   2.581  26.117 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.448582   0.794183   0.565    0.573    
# z.lag.1     -0.319587   0.049217  -6.493 4.91e-10 ***
#   tt          -0.002933   0.005678  -0.516    0.606    
# z.diff.lag   0.050771   0.064168   0.791    0.430    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.093 on 236 degrees of freedom
# Multiple R-squared:  0.1621,	Adjusted R-squared:  0.1514 
# F-statistic: 15.21 on 3 and 236 DF,  p-value: 4.394e-09
# 
# 
# Value of test-statistic is: -6.4934 14.1574 21.2186 
# 
# Critical values for test statistics: 
#   1pct  5pct 10pct
# tau3 -3.99 -3.43 -3.13
# phi2  6.22  4.75  4.07
# phi3  8.43  6.49  5.47



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

# 1b) run the ADF test allowing for an intercept + trend
adf_trend3 <- ur.df(
  y3,
  type       = "trend",    # intercept + time trend
  selectlags = "AIC"       # choose optimal lag length
)

# 1c) inspect the output
summary(adf_trend3)


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
# -45.161  -1.756   0.094   2.270  30.525 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.021780   0.841340   0.026    0.979    
# z.lag.1     -1.398209   0.093980 -14.878  < 2e-16 ***
#   tt          -0.002706   0.006094  -0.444    0.657    
# z.diff.lag   0.262852   0.063184   4.160 4.47e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.415 on 233 degrees of freedom
# Multiple R-squared:  0.5847,	Adjusted R-squared:  0.5793 
# F-statistic: 109.3 on 3 and 233 DF,  p-value: < 2.2e-16
# 
# 
# Value of test-statistic is: -14.8777 73.7838 110.6744 
# 
# Critical values for test statistics: 
#   1pct  5pct 10pct
# tau3 -3.99 -3.43 -3.13
# phi2  6.22  4.75  4.07
# phi3  8.43  6.49  5.47


#This time, we observe no more trend and drift, but the lag 1 coefficient is significant! 


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

# We observe that the new series does not have a drift and a trend anymore we then can try to fit an AR(1)


# End of Part 1
