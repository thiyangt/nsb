library(Mcomp)
library(forecast)
data(M3)

# Extract annual data
yearly_m3 <- subset(M3, "yearly")

# Simulate time series using automated ARIMA algorithm
set.seed(1)

# fit arima models to M3 data
arimaM3 <- lappy(yearly_m3, function(ts) auto.arima(as.ts(c(ts$x, ts$xx))))

simA <- lapply(arimaM3)
