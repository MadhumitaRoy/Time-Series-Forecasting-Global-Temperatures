library(tidyverse)
library(forecast)
library(fpp)

uk.met.data.file <- read.csv("UK Met Office_Data.csv",header = TRUE,sep = ",")
# Changing the column name to Temperature.diff. This is the relative temperature in Celcius
names(uk.met.data.file)[2] <- "Temperature.diff.C"
# Changing the temperature values to absolute temperatures
uk.met.data.file$Abs.temperature<- uk.met.data.file$Temperature.diff + 14
# Changing it into a time series with monthly frequency of 12
uk.met.data.file_ts <- ts(uk.met.data.file$Abs.temperature,start=1850,frequency = 12)

fit.stl <- stl(uk.met.data.file_ts, t.window=12, s.window="periodic") #decompose using STL (Season and trend using Loess)
plot(fit.stl)
plot(uk.met.data.file_ts)

# Looking at the data from 2000 to 2020
uk.met.data.file_ts_2000_2020 <- ts(uk.met.data.file$Abs.temperature,start=2000,end = 2020,frequency = 12)
plot(uk.met.data.file_ts_2000_2020)
#### ETS Model ####
# Create exponential smoothing models (ETS):

uk.met.data_AAZ <- ets(uk.met.data.file_ts, model="AAZ", damped=FALSE)
uk.met.data_MMZ <- ets(uk.met.data.file_ts, model="MMZ", damped=FALSE)

# Create their prediction "cones" for 969 months (till Dec 2100) into the future with quintile confidence intervals

uk.met.data_AAZ_pred <- forecast(uk.met.data_AAZ, h=969, level=c(0.8, 0.95))
uk.met.data_MMZ_pred <- forecast(uk.met.data_MMZ, h=969, level=c(0.8, 0.95))

# Compare the prediction "cones" visually
par(mfrow=c(1,2))
plot(uk.met.data_AAZ_pred, xlab="Year", ylab="Predicted Temperature in C")
plot(uk.met.data_MMZ_pred, xlab="Year", ylab="Predicted Temperature in C")

# Lets look at what our models actually are -- ETS
uk.met.data_AAZ
uk.met.data_MMZ

#### TBATS Model ####
#Create a trigonometric box-cox autoregressive trend seasonality (TBATS) model
uk.met.data_tbats <- tbats(uk.met.data.file_ts)
uk.met.data_tbats_pred <-forecast(uk.met.data_tbats, h=969, level=c(0.8, 0.95))
par(mfrow=c(1,1))
plot(uk.met.data_tbats_pred, xlab="Year", ylab="Predicted Temperature in C")

par(mfrow=c(1,3)) # Lets look at the three models with seasonality on one graph on the same scale
plot(uk.met.data_AAZ_pred, xlab="Year", ylab="Predicted Temperature in C", ylim=c(12,20))
plot(uk.met.data_MMZ_pred, xlab="Year", ylab="Predicted Temperature in C", ylim=c(12,20))
plot(uk.met.data_tbats_pred, xlab="Year", ylab="PPredicted Temperature in C", ylim=c(12,20))

# Lets look at what our models actually are -- TBATS
uk.met.data_tbats

## Calculating MAPE
f_MMM  <- function(y, h) forecast(ets(y, model="MMM"), h = h)
errors_MMM <- tsCV(uk.met.data.file_ts, f_MMM, h=25, window=1000)

## MAPE of MMM ETS model
mean(abs(errors_MMM/uk.met.data.file_ts), na.rm=TRUE)*100

f_AAA  <- function(y, h) forecast(ets(y, model="AAA"), h = h)
errors_MMM <- tsCV(uk.met.data.file_ts, f_MMM, h=25, window=1000)
## MAPE of MMM TBATS model
f_TBATS  <- function(y, h) forecast(tbats(y), h = h)
errors_TBATS <- tsCV(uk.met.data.file_ts, f_TBATS, h=25, window=1000)

mean(abs(errors_TBATS/uk.met.data.file_ts), na.rm=TRUE)*100

#Writing the csv file for ETS MMM model and TBATS model. TBATS is better than ETS
write.csv(uk.met.data_MMZ_pred, file = "UK MET Predicted temperatures MMM Model.csv")
write.csv(uk.met.data_tbats_pred, file = "UK MET Predicted temperatures TBATS Model.csv")

# Print the mean and confidence intervals for the MMZ model
uk.met.data_MMZ_pred

#### ARIMA MODEL ####
par(mfrow=c(1,1))
plot(uk.met.data.file_ts, xlab="Year", ylab="Uk temperatures")
plot(log(uk.met.data.file_ts), xlab="Year",ylab="log Uk temperatures")
plot(diff(log(uk.met.data.file_ts),12), xlab="Year",ylab="Annual change in monthly log A10 sales")

fit.stl <- stl(log(uk.met.data.file_ts), t.window=12, s.window="periodic")
plot(fit.stl)

# auto-correlation function
Acf(uk.met.data.file_ts,main="") # data "as is"
Acf(log(uk.met.data.file_ts),main="") # log-transformed data
Acf(diff(log(uk.met.data.file_ts),12),main="") # difference-12 log data
Pacf(diff(log(uk.met.data.file_ts),12),main="") 

fit.arima <- auto.arima(uk.met.data.file_ts,seasonal=TRUE)
fit.arima

par(mfrow=c(1,1))
plot(residuals(fit.arima))
Acf(residuals(fit.arima))

f_cast <- forecast(fit.arima,969,level=c(0.8,0.9, 0.95)) #969 stands for 969 months till Dec 2100
plot(f_cast,xlab = "Year",ylab = "Temperature(C)")

f_cast$mean<-f_cast$mean-14
f_cast$x<-f_cast$x-14
f_cast$lower<-f_cast$lower-14
f_cast$upper<-f_cast$upper-14

plot(f_cast)
write.csv(f_cast, file = "UK MET Predicted temperatures Arima Model.csv")

#Arima cross Validation and calculating MAPE
f_cast_arima  <- function(y,h)forecast(auto.arima(y),h=h)
errors_arima <- tsCV(uk.met.data.file_ts,f_cast_arima, h=25, window=1000)

mean(abs(errors_arima/uk.met.data.file_ts),na.rm=TRUE)*100

