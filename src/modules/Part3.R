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

################################### Question 8 ################################### 

theta1 <- coef(fit)["ma1"]
sigma2 <- fit$sigma2
# your point‐forecasts:
fc   <- forecast(fit, h = 2, level = 95)
mu1  <- fc$mean[1]   # Ŷ_{T+1}
mu2  <- fc$mean[2]   # Ŷ_{T+2}


Sigma <- sigma2 * matrix(c(
  1,        theta1,
  theta1,   1 + theta1^2
), nrow = 2, byrow = TRUE)

alpha <- 0.05
crit  <- qchisq(1 - alpha, df = 2)   # ≈ 5.991

install.packages("ellipse")   
library(ellipse)

el <- ellipse(Sigma, centre = c(mu1, mu2), level = 1 - alpha)
plot(el, type = "l", asp = 1,
     xlab = expression(X[T+1]), ylab = expression(X[T+2]),
     main = "95% Joint Confidence Ellipse")
points(mu1, mu2, pch = 19, col = "red")

