###• Part II : ARMA models

################################### Question 4 ################################### 
#We look at the ACF/PACF 


# suppose your cleaned series is in data_ts$differenced 
x <- na.omit(data_ts$differenced)  

# show both plots side by side
ggAcf(x, lag.max=24) + ggtitle("ACF of X_t")
ggPacf(x, lag.max=24) + ggtitle("PACF of X_t")


#ACF cuts off after lag q ⇒ suggests MA(q). Then we can try MA(2)

#PACF cuts off after lag p ⇒ suggests AR(p) Then we can try AR(4)



# 1 fit MA(2) with no constant (we removed trend/drift)
fit_ma2 <- Arima(
  x,
  order            = c(0, 0, 2),
  include.constant = FALSE
)

summary(fit_2)

 
# Series: x 
# ARIMA(0,0,2) with zero mean 
# 
# Coefficients:
#   ma1      ma2
# -0.2731  -0.3575
# s.e.   0.0589   0.0582
# 
# sigma^2 = 37.98:  log likelihood = -772.99
# AIC=1551.97   AICc=1552.07   BIC=1562.4
# 
# Training set error measures:
#   ME     RMSE      MAE     MPE     MAPE      MASE       ACF1
# Training set -0.5408261 6.136934 3.280441 104.347 218.9111 0.5953705 0.02815581

coefs <- coef(fit_ma2)
ses   <- sqrt(diag(fit_ma2$var.coef))
zvals <- coefs / ses
pvals <- 2*(1 - pnorm(abs(zvals)))
round(data.frame(Estimate=coefs, Std.Error=ses,  z=zvals, p=pvals), 4)

#Estimate Std.Error       z p
#ma1  -0.2731    0.0589 -4.6382 0
#ma2  -0.3575    0.0582 -6.1417 0


#We computed the z-statistics by dividing each estimate by its standard error. Both ∣z∣>>1.96, so ma1 and ma2 are highly significant at the 1% level.

#Now the invertibility condition (the modules of the roots of the MA polynom must be outside the unit circle):

# 1) Extract the MA coefficients

theta1 <- coef(fit_ma2)["ma1"]
theta2 <- coef(fit_ma2)["ma2"]

# 2) Form the polynomial 1 + theta1*z + theta2*z^2

roots <- polyroot(c(1, theta1, theta2))

# 3) Check their modules
abs_roots <- Mod(roots)
abs_roots

#Finally we have to check is the residuals act as White Noise: 



checkresiduals(fit_ma2)








# 2 fit AR(4) with no constant (we removed trend/drift)
fit_ar4 <- Arima(
  x,
  order            = c(4, 0, 0),
  include.constant = FALSE
)


summary(fit_ar4)


# Series: x 
# ARIMA(4,0,0) with zero mean 
# 
# Coefficients:
#   ar1      ar2      ar3      ar4
# -0.2004  -0.3250  -0.1894  -0.1559
# s.e.   0.0639   0.0639   0.0637   0.0635
# 
# sigma^2 = 38.92:  log likelihood = -774.82
# AIC=1559.64   AICc=1559.9   BIC=1577.02
# 
# Training set error measures:
#   ME     RMSE      MAE      MPE     MAPE      MASE        ACF1
# Training set -0.3711186 6.186071 3.317631 101.1735 178.4067 0.6021201 -0.02231636

coefs <- coef(fit_ar4)
ses   <- sqrt(diag(fit_ar4$var.coef))
zvals <- coefs / ses
pvals <- 2*(1 - pnorm(abs(zvals)))
round(data.frame(Estimate=coefs, Std.Error=ses,  z=zvals, p=pvals), 4)

# Estimate Std.Error       z      p
# ar1  -0.2004    0.0639 -3.1373 0.0017
# ar2  -0.3250    0.0639 -5.0851 0.0000
# ar3  -0.1894    0.0637 -2.9722 0.0030
# ar4  -0.1559    0.0635 -2.4556 0.0141


#Again all the coeff are highly significant.

#Now the invertibility condition (the modules of the roots of the MA polynom must be outside the unit circle):


# 1) Extract the AR coefficients

phi1 <- coef(fit_ar4)["ar1"]
phi2 <- coef(fit_ar4)["ar2"]
phi3 <- coef(fit_ar4)["ar3"]
phi4 <- coef(fit_ar4)["ar4"]

# 2) Form the polynomial 1 + phi1*z + phi2*z^2 + ... 

roots <- polyroot(c(1, phi1, phi2, phi3, phi4))

# 3) Check their modules
abs_roots <- Mod(roots)
abs_roots

#Check the residuals 

checkresiduals(fit_ar4)

#Precision the Ljung-Box test : 

#H₀: “no autocorrelation up to lag m.”

#If p-value > 0.05, you fail to reject ⇒ residuals behave like white noise.

################################### Question 5 ################################### 

# Since I have differenced the series once to obtain “differenced,” an ARIMA(0,1,2) or ARIMA(4,1,0) model could fit.

## AR(4) AIC=1559.64   AICc=1559.9   BIC=1577.02
## MA(2) AIC=1551.97   AICc=1552.07   BIC=1562.4

#Then we can choose ARIMA(0,1,2)



# END part 2