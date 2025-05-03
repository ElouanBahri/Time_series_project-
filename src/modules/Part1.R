# PART 1: Exploratory analysis and stationarity transformation for the Building Construction IPI series
# ==============================================================================
# 1. Description of the series:
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

library(zoo)

data_ts %>%
  mutate(
    roll_mean  = zoo::rollmean(Index, 24, fill = NA),
    roll_sd    = zoo::rollapply(Index, 24, sd, fill = NA)
  ) %>%
  ggplot(aes(x = ym)) +
  geom_line(aes(y = Index)) +
  geom_line(aes(y = roll_mean),   color = "blue") +
  geom_line(aes(y = roll_sd * 10), color = "red") 

# Even if the augmented DCF says that the series is stationary, it is evident that there is a dertministic trend downward we can test it :

# add a simple time counter 1,2,…,T
data_trend <- data_ts %>%
  mutate(t = row_number())

# regress index on t
fit_lin <- lm(Index ~ t, data = data_trend)
summary(fit_lin)

#The p-value of t is very low then t is statistically significant and there is well a dertiministic trend downward


#Correction of the deterministic trend 

data_diff <- data_ts %>%
  mutate(d1 = Index - lag(Index))

# Test for unit root on the series
adf_log <- adf.test(na.omit(data_diff$d1), alternative = "stationary")
print(adf_log)  # If p-value > 0.05, The series is non-stationary





## I am here 





# 7. Test for unit root on log series
adf_log <- adf.test(data_ts$Index, alternative = "stationary")
print(adf_log)  # If p-value > 0.05, log series is non-stationary

# 8. Difference the log series if non-stationary
#    (first difference approximates monthly growth rates)
data_ts <- data_ts %>%
  mutate(dlog_index = difference(log_index))

# 9. Test for stationarity on differenced log series
adf_diff <- adf.test(na.omit(data_ts$dlog_index), alternative = "stationary")
print(adf_diff)  # Expect p-value < 0.05 indicating stationarity

# 10. Plot log series and differenced log series
p1 <- ggplot(data_ts, aes(x = month, y = log_index)) +
  geom_line() +
  labs(title = "Log-transformed IPI - Building Construction",
       x = "Month", y = "Log(Index)") +
  theme_minimal()

p2 <- ggplot(data_ts, aes(x = month, y = dlog_index)) +
  geom_line() +
  labs(title = "First Difference of Log IPI (Monthly Growth Rate)",
       x = "Month", y = "Delta Log(Index)") +
  theme_minimal()

print(p1)
print(p2)

# End of Part 1
