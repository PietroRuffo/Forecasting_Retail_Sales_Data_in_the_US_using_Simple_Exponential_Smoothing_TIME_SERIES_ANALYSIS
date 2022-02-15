#######################################################
# This file loads monthly sales per day.              #
# This data has been adjusted for inflation.          #
# We will use data to forecast the next two years     #
# of sales per day.                                   #
#                                                     #
# Created by: Pietro Ruffo 27/11/21                   #
#######################################################

# Clear all variables in workspace
rm(list=ls())
######################################
#   P R E   P R O C E S S I N G      #
######################################

# Load libraries
library(tseries)
library(ggplot2)
library(ggthemes)
library(forecast)
library(zoo)
library(fpp2)

# Set working directory and Load the data
setwd("C:/Users/Utente/Desktop/UNI/ERASMUS LISBON/Processos de Previsão e Decisão/Data_Project_Pietro_Ruffo")
data<- read.csv("real_sales_per_day.csv")

attach(data)#to make the code more compact

# Check for missing values
sum(is.na(data)) #no missings

head(data) # the time series begins on the first of January of 1992
tail(data) # there are 312 observations, up to the first of December 2017 

# Declare this as time series data
Y<- ts(data[,2], start=c(1992,1),frequency=12) # monthly observations

#########################################################
#  E X P L O R A T O R Y   D A T A    A N A L Y S I S   #
#########################################################

# Observe the graphical representation of the series
autoplot(Y) + 
  ggtitle("Real US Retail Sales per Day") + 
  ylab("Millions of 2017 Dollars")

# The graph shows a seasonal pattern evolving on a positive trend. 

# TREND AND SEASONALITY DETECTION WITH OLS METHOD  #

# SEASONALITY DETECTION (seasonal dummies)

# Let's study the seasonal component through OLS method
# let's suppose the time series is defined by an additive model containing 
# the seasonal component plus the erratic term

# let's reorganize data so that it is possible to clearly identify the months
yy<-rep(1:12,26) #26 full years in total

data_seasonality<- cbind(yy,Y)

#creating the dummies for each month
d1<-ifelse(yy==1,1,0) #if yy=1 (January) then d1=1, else d1=0
d2<-ifelse(yy==2,1,0) #February
d3<-ifelse(yy==3,1,0) #March
d4<-ifelse(yy==4,1,0) #April
d5<-ifelse(yy==5,1,0) #May
d6<-ifelse(yy==6,1,0) #June
d7<-ifelse(yy==7,1,0) #July
d8<-ifelse(yy==8,1,0) #August
d9<-ifelse(yy==9,1,0) #September
d10<-ifelse(yy==10,1,0) #October
d11<-ifelse(yy==11,1,0) #November
d12<-ifelse(yy==12,1,0) #December

# The resulting matrix is the following
D<-cbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12);D
head(D)
tail(D)
data_seasonality1<-cbind(yy,Y,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12);data_seasonality1 

seas_data<-lm(Y~d1+d2+d3+d4+d5+d6+d7+d8+d9+d10+d11+d12-1) #Intercept deleted not to fall in the dummy trap
summary(seas_data)

# Data seasonally adjusted
res_des<-seas_data$residuals 
min(res_des)
res_des<-seas_data$residuals + abs(min(res_des))

# transformig the residuals in a time series object to be plugged into the autoplot function
res_desT<- ts(res_des, start=c(1992,1),frequency=12) 

autoplot(res_desT) + 
  ggtitle("Seasonally adjusted Series") + 
  ylab("Millions of 2017 Dollars")

ggseasonplot(res_desT) + 
  ggtitle("Seasonal Adjusted Plot: Real US Retail Sales per Day") + 
  ylab("Millions of 2017 Dollars")

# TREND DETECTION

# Adjusted Rsquared criterion
t<-seq(1, length(Y),1)
poly1<-lm(Y~poly(t,1,raw=T)); summary(poly1) # Adjusted R-squared:  0.6532
poly2<-lm(Y~poly(t,2,raw=T)); summary(poly2) # Adjusted R-squared:  0.6905
poly3<-lm(Y~poly(t,3,raw=T)); summary(poly3) # Adjusted R-squared:  0.7313 
poly4<-lm(Y~poly(t,4,raw=T)); summary(poly4) # Adjusted R-squared:  0.7446 <- 4th degree polynomial 
poly5<-lm(Y~poly(t,5,raw=T)); summary(poly5) # Adjusted R-squared:  0.7441 

#  fitting a 4th degree polynomial
fit_trend<-fitted(poly4)
ts.plot(Y, fit_trend, gpars = list(col=c(1,6),main="Observed vs Fitted values") )

# Series trend adjusted
res_det<-poly4$residuals

res_detT<- ts(res_det, start=c(1992,1),frequency=12) 

autoplot(res_detT) + 
  ggtitle("Trend adjusted Series") + 
  ylab("Millions of 2017 Dollars")

#season plot for trend adjusted series
ggseasonplot(res_detT) + 
  ggtitle("Seasonal Plot: Real US Retail Sales per Day") + 
  ylab("Millions of 2017 Dollars")
#Let's look at another seasonal plot, the subseries plot
ggsubseriesplot(res_detT)

# Trend plus Seasonality detection
tr_ss_Y<-lm(Y ~ d1+d2+d3+d4+d5+d6+d7+d8+d9+d10+d11+d12-1 + poly(t,4,raw=T)); summary(tr_ss_Y)
fit_ss_tr<-fitted(tr_ss_Y)
ts.plot(Y, fit_ss_tr, gpars = list(col=c(1,6),main="Observed vs Fitted values") )

# Fitted value
fit_tr_ss_Y<-fitted(tr_ss_Y); fit_tr_ss_Y
# Plotting fitted and observed values
ts.plot(Y, fit_tr_ss_Y, gpars = list(col=c(1,6),main="Observed vs Fitted values") )

##########################
# F O R E C A S T I N G  #
##########################

#ACF dacays over time
ggAcf(Y, lag=60)
ggAcf(res_des, lag=60)

# Gadner-McKenzie METHOD

# Let's forecast using the series filtred of seasonality
# plugging different value for beta
fc1<- holt(res_des,damped = T, phi=NULL, h=24);summary(fc1) #alpha=0.449, beta=0.0315,MAE=180.0497
#AIC=5178.622 and BIC=5201.080 
fc2<-holt(res_des, damped = T, phi=NULL, beta=0.2, h=24);summary(fc2) #MAE=186.2401 
fc3<-holt(res_des, damped = T, phi=NULL, beta=0.4, h=24);summary(fc3) #MAE= 192.8345
fc4<-holt(res_des, damped = T, phi=NULL, beta=0.8, h=24);summary(fc4) #MAE= 269.2096

fc1$fitted
#Observed vs fitted values Gadner_McKenzie Method
plot(fc1$fitted)
ts.plot(res_des, fc1$fitted, gpars = list(col=c(1,6),main="Gadner-McKenzie  Method: Observed vs Fitted Values ") )

#cost function beta parameter
beta_holt <- cbind(0.0315,0.2,0.4,0.8)
MAE_holt<-cbind(178.8047,186.2401,192.8345,269.2096)

plot(beta_holt,MAE_holt ,type = "b", col = 2 , lwd = 3, pch = 1, main = "Cost Function Beta parameter Gadner-McKenzie")
#the closer to 0 the smaller the MAE
autoplot(fc1) + 
  ggtitle("Gadner-McKenzie Forecasting") + 
  ylab("Millions of 2017 Dollars")

#HOLT LINEAR VS DAMPED
fc <- holt(res_desT, h=24);summary(fc)
fc2 <- holt(res_desT, damped=TRUE, phi = 0.9457, h=24)
autoplot(res_desT) +
  autolayer(fc, series="Holt's method", PI=FALSE) +
  autolayer(fc2, series="Gadner-McKenzie's method", PI=FALSE) +
  ggtitle("Real US Retail Sales per Day") + 
  ylab("Millions of 2017 Dollars")
guides(colour=guide_legend(title="Forecast"))

# Let's forecast  the initial data Y containing both a trend and an additive seasonality

# HOLT - WINTERS
#testing for different values of alpha
holt_win<-hw(Y,h=24, seasonal = "additive");summary(holt_win)#MAE= 171.0885
holt_win2<-hw(Y,h=24, seasonal = "additive", alpha = 0.2);summary(holt_win2)# MAE=175.4056
holt_win3<-hw(Y,h=24, seasonal = "additive", alpha = 0.6);summary(holt_win3)# MAE = 175.1954
holt_win4<-hw(Y,h=24, seasonal = "additive", alpha = 0.8);summary(holt_win4)# MAE= 189.7376

#cost function beta parameter
alpha_holtw <- cbind(0.2,0.402,0.6,0.8)
MAE_holtw<-cbind(175.4056,171.0885,175.4056,189.7376)

plot(alpha_holtw,MAE_holtw ,type = "b", col = 2 , lwd = 3, pch = 1, main = "Cost Function Alpha parameter Holt-Winters")
#the closer to 0 the smaller the MAE

#Observed vs fitted values Holt Winters
fitted<-holt_win$fitted
plot(fitted)
ts.plot(Y, fitted, gpars = list(col=c(1,6),main="Holt-Winters  Method: Observed vs Fitted Values ") )

#Holt-Winters FORECASTING
hw_Y_forecasting <- forecast(holt_win, h=24)
hw_Y_forecasting
autoplot(hw_Y_forecasting) + 
  ggtitle("Holt-Winters Method") + 
  ylab("Millions of 2017 Dollars")

####################
####EXTRA IDEAS#####
####################

#Let's look at another seasonal plot, the subseries plot
ggsubseriesplot(res_desT)

#alternative method to remove seasonality
decomp <- decompose(Y)
autoplot(decomp)
res_des1 <- Y  - decomp$seasonal
ts.plot(res_des1)
# trend detection by differences method (graphical method)
d1_Y<-diff(Y, lag=1, differences = 1)
autoplot(d1_Y) + 
  ggtitle("First Order Difference") + 
  ylab("Millions of 2017 Dollars")


d2_Y<-diff(Y, lag=1, differences = 2)
autoplot(d2_Y) + 
  ggtitle("Second Order Difference") + 
  ylab("Millions of 2017 Dollars")

d3_Y<-diff(Y, lag=1, differences = 3)
autoplot(d3_Y) + 
  ggtitle("Third Order Difference") + 
  ylab("Millions of 2017 Dollars")

d4_Y<-diff(Y, lag=1, differences = 4)
autoplot(d4_Y) + 
  ggtitle("Fourth Order Difference") + 
  ylab("Millions of 2017 Dollars")

# It might be a fourth order polynomial
# Series flitred of Trend and Seasonality
filtred_Y<- tr_ss_Y$residuals
filtred_Y<- ts(filtred_Y, start=c(1992,1),frequency=12) 

autoplot(filtred_Y) + 
  ggtitle("Series adjusted for Trand and Seasonality") + 
  ylab("Millions of 2017 Dollars")

#season plot for Series flitred of Trend and Seasonality
ggseasonplot(filtred_Y) + 
  ggtitle("Seasonal Plot: Real US Retail Sales per Day") + 
  ylab("Millions of 2017 Dollars")
#Let's look at another seasonal plot, the subseries plot
ggsubseriesplot(filtred_Y)

#SIMPLE EXPONENTIAL SMOOTHING

#using the the series filtred of both trend and seasonality
#let's try different value of alpha
ses1<- ses(filtred_Y, h=24, alpha = 0.425);summary(ses1) #MAE=169.584
ses2<- ses(filtred_Y, h=24, alpha = 0.2);summary(ses2) #MAE=178.6794 
ses3<- ses(filtred_Y, h=24, alpha = 0.6);summary(ses3) #MAE= 175.7273
ses4<- ses(filtred_Y, h=24, alpha = 0.8);summary(ses4) #MAE= 190.1511

alpha_ses <- cbind(0.2,0.425,0.6,0.8)
MAE<-cbind(178.6794,169.584,175.7273,190.1511)

plot(alpha_ses,MAE ,type = "b", col = 2 , lwd = 3, pch = 1, main = "Cost Function Alpha parameter SES")
lines(alpha_ses,MAE)

#SES forecasting
ses_forecasting <- forecast(ses1, h=24)
ses_forecasting
autoplot(ses_forecasting) + 
  ggtitle("SES Forecasting") + 
  ylab("Millions of 2017 Dollars")

# HOLT METHOD

# Let's forecast using the series filtred of seasonality
# plugging different value for beta
hfc1<- holt(res_des,damped = F, h=24);summary(hfc1) #alpha=0.449, beta=1e-04 ,MAE=178.8047
#AIC= 5176.671 , BIC=5195.386 , MAE = 178.8047
hfc2<-holt(res_des, damped = F, beta=0.2, h=24);summary(hfc2) #MAE=195.7573 
hfc3<-holt(res_des, damped = F, beta=0.4, h=24);summary(hfc3) #MAE= 202.6984
hfc4<-holt(res_des, damped = F, beta=0.8, h=24);summary(hfc4) #MAE= 288.2785

#Observed vs fitted values Holt Method
plot(fc1$fitted)
ts.plot(res_des, hfc1$fitted, gpars = list(col=c(1,6),main="Gadner-McKenzie  Method: Observed vs Fitted Values ") )

#cost function beta parameter
hbeta_holt <- cbind(1e-04,0.2,0.4,0.8)
hMAE_holt<-cbind(178.8047,195.7573,202.6984,288.2785)

plot(hbeta_holt,hMAE_holt ,type = "b", col = 2 , lwd = 3, pch = 1, main = "Cost Function Beta parameter HOLT")
#the closer to 0 the smaller the MAE
autoplot(hfc1) + 
  ggtitle("Holt Forecasting") + 
  ylab("Millions of 2017 Dollars")
