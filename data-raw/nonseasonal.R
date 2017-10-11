library(Mcomp)
library(forecast)
data(M3)

# Extract annual data
yearly_m3 <- subset(M3, "yearly")

# Simulate time series using automated ARIMA algorithm
set.seed(1)

# fit arima models to M3 data
arimaM3 <- lapply(yearly_m3, function(ts) auto.arima(as.ts(c(ts$x, ts$xx))))

# simulated 1000 series from the models in arimaM3.
# simA is a list, and each element in the list is a matrix of 1000 columns and
# number of raws in the matrix equals to the full length(training+test) of the corresponding series in M3
# simA returns a list of matrices each contains 1000 columns. 645 matrices are in the list.
simA <- lapply(arimaM3,
												function(temp){
												length_series <- length(temp)
												mat <- ts(matrix(0, ncol=1000, nrow=length_series))
												for(i in 1:1000)
												mat[,i] <- simulate(temp, nsim=length_series)
												return (mat)})

# Finding the best forecasting model for each simulated series

tail.ts <- function(data, n){
	data <- as.ts(data)
	window(data, start=tsp(data)[2]-(n-1)/frequency(data))
}

head.ts <- function(data, n){
	data <- as.ts(data)
	window(data, end=n)
}

best_method <- function(x, xx){ # x - time series, xx-number of points for the test set

	training_len <- length(x) - xx
	#temp <- as.ts(temp)
	testing <- tail.ts(x, xx)
	training <- head.ts(x, training_len)
	h <- length(testing)

	#auto.arima
	fit.arima <- auto.arima(training,max.p=3, max.q=3)
	forecast.arima <- forecast(fit.arima,h)
	autoArima_mase <- accuracy(forecast.arima,testing)["Test set","MASE"]

	#fit random walk
	fit.RW <- rwf(training,drift=FALSE)
	RW_mase <- accuracy(fit.RW,testing)["Test set","MASE"]

	#fit random walk with drift
	fit.RWD <- rwf(training,drift=TRUE)
	RWD_mase <- accuracy(fit.RWD,testing)["Test set","MASE"]

	#fit white noise
	fit.WN <- Arima(training,order=c(0,0,0))
	forecast.WN <- forecast(fit.WN,h)
	WN_mase <- accuracy(forecast.WN,testing)["Test set","MASE"]

	#fit ets
	fit.ets <- ets(training)
	forecast.ets <- forecast(fit.ets,h)
	ets_mase <- accuracy(forecast.ets,testing)["Test set","MASE"]

	#list_fits <- list(fit.arima, fit.RW, fit.RWD, fit.WN, fit.ets)
	mase_val <- c(autoArima_mase, RW_mase, RWD_mase, WN_mase, ets_mase)
	min_mase <- which.min(mase_val)

	if (min_mase==1) {
		best_method <-  forecast.arima$method
	} else if (min_mase==2) {
		best_method <- fit.RW$method
	} else if (min_mase==3) {
		best_method <-  fit.RWD$method
	} else if (min_mase==4) {
		best_method <-  forecast.WN$method
	} else {
		best_method <- forecast.ets$method
	}
return(best_method)
}

BESTsimA <- lapply(simA, apply, 2, best_method, xx=6)

# Simulate series with ETS function

#etsM3 <- lapply(yearly_m3, function(ts) ets(as.ts(c(ts$x, ts$xx))))

# simulated 1000 series from the models in arimaM3.
# simA is a list, and each element in the list is a matrix of 1000 columns and
# number of raws in the matrix equals to the full length(training+test) of the corresponding series in M3
# simA returns a list of matrices each contains 1000 columns. 645 matrices are in the list.
set.seed(2)
simE <- lapply(etsM3,
							 function(temp){
							 	length_series <- length(temp)
							 	matE <- ts(matrix(0, ncol=1000, nrow=length_series))
							 	for(i in 1:1000)
							 		matE[,i] <- simulate(temp, nsim=length_series)
							 	return (matE)})

# simulate function note:
#Future=T, the simulated observations are conditional on the
#historical observations. In other words, they are possible future sample paths
#of the time series. But if future=FALSE, the historical data are ignored, and the
#simulations are possible realizations of the time series models that are
#not connected to the original data.
#

BESTsimE <- lapply(simE, apply, 2, best_method, xx=6)

# Things that should move to the package
#simA- simulated time series by using auto.arima
#simE - simulated time series using ets
#BESTsimA - best forecasting method for each simulated time series in simA
#BESTsimE - best forecasting method for each simulated time series in simE

devtools::use_data(simA)
devtools::use_data(BESTsimA)

devtools::use_data(simE)
devtools::use_data(BESTsimE)

##########################################################
#
# Best forecasting method for each yearly M3 series and MASE calculated from each method
#(auto.arima, ets, rwd, rw, wn)
#
###########################################################

bestFORm3Annual <- lapply(yearly_m3,
													function(temp){
														training <- temp$x
														testing <- temp$xx
														h <- length(testing)

														# auto.arima
														fitArima <- auto.arima(training)
														forecastArima <- forecast(fitArima, h)
														ARIMA <- accuracy(forecastArima,testing)[2,6]


														# ETS models
														fitEts <- ets(training)
														forecastEts <- forecast(fitEts, h)
														ETS <- accuracy(forecastEts,testing)[2,6]

														# WN

														fitWN <- Arima(training,order=c(0,0,0))
														forecastWN <- forecast(fitWN, h)
														WN <- accuracy(forecastWN,testing)[2,6]

														# RW

														fitRW <- Arima(training,order=c(0,1,0))
														forecastRW <- forecast(fitRW, h)
														RW <- accuracy(forecastRW,testing)[2,6]

														# RD
														fitRWD <- Arima(training,order=c(0,1,0), include.drift=T)
														forecastRWD <- forecast(fitRWD,h)
														RWD <- accuracy(forecastRWD,testing)[2,6]


														MaseM3 <- data.frame(ARIMA, ETS, WN, RW, RWD)

														names(MaseM3) <- c("auto.arima", "ets", "WN", "RW", "RWD")
														MaseM3_min <- which.min(MaseM3)

														if (MaseM3_min==1) {
															MaseM3$bestMethod <- as.character(fitArima)
														} else if (MaseM3_min==2) {
															MaseM3$bestMethod <- as.character(fitEts)
														} else if (MaseM3_min==3) {
															MaseM3$bestMethod <- as.character(fitWN)
														}  else if (MaseM3_min==4) {
															MaseM3$bestMethod <- as.character(fitRW)
														}else
															MaseM3$bestMethod <- as.character(fitRWD)

														return(MaseM3)


													})
maseM3df <- do.call("rbind", bestFORm3Annual) # Combine all dataframes into one
library(dplyr)
maseM3df <- dplyr::add_rownames(maseM3df, "ID")

devtools::use_data(maseM3df)

# Finding the best MASE from each method (auto.arima, ets, rw, rwd,wn)for each time series in
# the M1 competition
data(M1)

# Extract annual data
yearly_m1 <- subset(M1, "yearly")

bestFORm1Annual <- lapply(yearly_m1,
													function(temp){
														training <- temp$x
														testing <- temp$xx
														h <- length(testing)

														# auto.arima
														fitArima <- auto.arima(training)
														forecastArima <- forecast(fitArima, h)
														ARIMA <- accuracy(forecastArima,testing)[2,6]


														# ETS models
														fitEts <- ets(training)
														forecastEts <- forecast(fitEts, h)
														ETS <- accuracy(forecastEts,testing)[2,6]

														# WN

														fitWN <- Arima(training,order=c(0,0,0))
														forecastWN <- forecast(fitWN, h)
														WN <- accuracy(forecastWN,testing)[2,6]

														# RW

														fitRW <- Arima(training,order=c(0,1,0))
														forecastRW <- forecast(fitRW, h)
														RW <- accuracy(forecastRW,testing)[2,6]

														# RD
														fitRWD <- Arima(training,order=c(0,1,0), include.drift=T)
														forecastRWD <- forecast(fitRWD,h)
														RWD <- accuracy(forecastRWD,testing)[2,6]


														MaseM3 <- data.frame(ARIMA, ETS, WN, RW, RWD)

														names(MaseM3) <- c("auto.arima", "ets", "WN", "RW", "RWD")
														MaseM3_min <- which.min(MaseM3)

														if (MaseM3_min==1) {
															MaseM3$bestMethod <- as.character(fitArima)
														} else if (MaseM3_min==2) {
															MaseM3$bestMethod <- as.character(fitEts)
														} else if (MaseM3_min==3) {
															MaseM3$bestMethod <- as.character(fitWN)
														}  else if (MaseM3_min==4) {
															MaseM3$bestMethod <- as.character(fitRW)
														}else
															MaseM3$bestMethod <- as.character(fitRWD)

														return(MaseM3)


													})
maseM1df <- do.call("rbind", bestFORm1Annual) # Combine all dataframes into one
library(dplyr)
maseM1df <- dplyr::add_rownames(maseM1df, "ID")
devtools::use_data(maseM1df)
