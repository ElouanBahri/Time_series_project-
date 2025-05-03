###• Part III : Prediction

#Some prediction 

library(forecast)



fit <- Arima(y, order = c(0,1,2), include.constant=FALSE)

# 1) Forecast the next H points, e.g. next 12 months
H <- 12
fc  <- forecast(fit, h = H, level = c(80, 95)) 


# 2) Inspect the result
print(fc)
#   - fc$mean: the point forecasts for t = T+1 … T+H
#   - fc$lower, fc$upper: the 80% and 95% interval bounds

# 3) Plot the forecasts against your historical series
autoplot(fc) +
  autolayer(fitted(fit), series = "Fitted") +
  labs(
    title = "Forecasts from ARIMA Model",
    x     = "Time",
    y     = "Series"
  ) +
  theme_minimal()
