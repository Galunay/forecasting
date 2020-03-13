#Preambula
library(tseries)
library(forecast)
library(ggplot2)
library(seasonal)
library(tidyverse)
library(rio)
library(xts)
library(gridExtra)
library(corrplot)
library(lmtest)
library(vars)
library(BigVAR)
library(tsDyn)
library(urca)

setwd("~/Documents/GitHub/bank_inst_has/metrics")
#setwd('C:\\Users\\123\\Desktop\\metricsGalya\\met\\')

#Data import
dd=import('data_for_HW1.xlsx')
glimpse(dd)
colSums(is.na(dd))
date=dplyr::select(dd, Time)
dd=dplyr::select(dd, -Time)

var.data  <- ts(dd[,1:5], start = c(2010, 1,1), frequency = 12)
#var.data[,1] <- log(var.data[,1])

#Data Visualisation, Stationarity 
plot(var.data
  )
##all ts are non stationary in row data (нужно сделать фасетку с acf pacf, попробовать gg)
for (i in 1:5){
  a<-ggAcf(var.data[,i])
  b<-ggPacf(var.data[,i])
  print(a)
  print(b)
}

#ADF testing
##all of them are non stationary
for (i in 1:5){
  c<-adf.test(var.data[,i])
  print(c)
}

#How should we modify process? Box Cox procedure
#Read here 
#https://www.statisticshowto.datasciencecentral.com/box-cox-transformation/

ggplot(dd, aes(x=log(var.data[,1])))+geom_histogram()

for (i in 1:5){
  d<-BoxCox.lambda(var.data[,i])
  print(d)
}


#Data transformation for ARIMA only
#Actually autoarima can work with non stat and work successfully as diff the data by itself
var.data[,1] <- log(var.data[,1])
var.data[,2] <- log(var.data[,2])
var.data[,3] <- (sqrt(var.data[,3]))^(-1)
var.data[,4] <- log(var.data[,4])
plot(var.data)

#Split into train and test data
train.var <- window(var.data, start=c(2010,1), 
                    end=c(2016,12), frequency=12) #The first train sample
test.var  <- window(var.data, start=c(2017,1), frequency=12) # whole test sample
hor  <- 12
Nsmp <- nrow(test.var)

# train_log=log(train)
# plot(train_log)
# ggAcf(train_log)
# adf.test(diff(train_log))

#Structure breaks at least I've tried to do it
##Иван Павлович, тут Галю Пажитнову понесло куда-то в дебри, хз что она тут творит
# library(tsDyn)
# dts<-setar(diff(train),1)
# plot(dts)

require(strucchange)
for (i in 1:5){
  bp.train <- breakpoints(var.data[,i] ~ 1)
  summary(bp.train)
  plot(bp.train)
  plot(var.data[,i])
  lines(bp.train)
  ## confidence intervals
  ci.train <- confint(bp.train)
  ci.train
  lines(ci.train)
}
# нужно дамми на конец 2014 - начало 2015 год, там была смена режима курса 
# и переход на инфляционное таргетирование, 
# правила можно найти на сайте цб
dummy <- var.data*0 
window(dummy, start = c(2014,12)) <- 1

#Some useful plots
ggseasonplot(train.var[,1]) #Seasonal subseries plot
ggseasonplot(train.var[,1], polar = T) #Seasonal polar plot

#Cointegration VECM prep
summary(ca.jo(var.data))
p1<-predict(vec2var(ca.jo(var.data), r=1), h=12)
df <- ts(p1$fcst$Deposits_ind[,1], start = c(2020,1))
autoplot(var.data[,1], series = 'Actual Data') + autolayer(df, series = 'Forecast')

##Naïve forecasting
plot(train.var)
n1<-naive(train.var[,1], h=12)
autoplot(n1)

##Models
models <- c('arima', 'ets', 'arimax', 'seas', 'var', 'var_lasso','vecm', 'rwd')
e <- f <-fL<-fH <- lapply(1:length(models), 
            function(i) {matrix(NA, hor, (Nsmp - hor))})
names(e) <- names(f) <-names(fL) <-names(fH) <- models

# NB! Initialize TT before loop!!!
TT   <- nrow(train.var)
nn = ncol(var.data)
criterion <- 1
for (i in 1:(Nsmp - hor)){
  cat(paste0('Sample ', TT, '\n'))
  # Form train samples
  train.yi <- var.data[1:TT, 1]
  train.xi <- var.data[1:TT, 2:5]
  dtrain.i <- dummy[1:TT]
  
  # Form test samples
  test.yi <- var.data[(TT+1):(TT+hor), 1]
  test.xi <- var.data[(TT+1):(TT+hor), 2:5]
  dtest.i <- dummy[(TT+1):(TT+hor)]
  
  # ARIMA
  m.arima <- auto.arima(train.yi, xreg = dtrain.i)
  f[['arima']][, i] <- forecast( m.arima , xreg = dtest.i, h=hor)$mean
  fL[['arima']][, i] <- forecast( m.arima , xreg = dtest.i, h=hor)$lower[,2]
  fH[['arima']][, i] <- forecast( m.arima , xreg = dtest.i, h=hor)$upper[,2]
  
  # ARIMAX
  xtr <- cbind(dtrain.i, train.xi); colnames(xtr) <- NULL
  xte <- cbind(dtest.i, test.xi); colnames(xte) <- NULL
  
  m.arimax <- auto.arima(train.yi, xreg = xtr)
  f[['arimax']][, i] <- forecast( m.arimax ,
                                  xreg = xte,
                                  h=hor)$mean
  fL[['arimax']][, i] <- forecast( m.arimax ,
                                   xreg = xte,
                                   h=hor)$lower[,2]
  fH[['arimax']][, i] <- forecast( m.arimax ,
                                   xreg = xte,
                                   h=hor)$upper[,2]
  
  # X13 SEATS
  seasD <- seas(ts(train.yi, frequency = 12, start = c(2010, 1)))
  seasComponent = seasD$series$s10
  trendComponent = seasD$series$s12
  
  f[['seas']][, i] <- snaive(seasComponent, h = hor)$mean *  
    rwf(trendComponent, drift = TRUE,  h = hor)$mean
  fL[['seas']][, i] <- snaive(seasComponent, h = hor)$mean *  
    rwf(trendComponent, drift = TRUE, h = hor)$lower[,2]
  fH[['seas']][, i] <- snaive(seasComponent, h = hor)$mean *  
    rwf(trendComponent, drift = TRUE,  h = hor)$upper[,2]
  
  # ETS
  m.ets <- ets(train.yi)
  f[['ets']][, i] <- forecast(m.ets, h=hor)$mean
  fL[['ets']][, i] <- forecast( m.ets,
                                   h=hor)$lower[,2]
  fH[['ets']][, i] <- forecast( m.ets,
                                   h=hor)$upper[,2]
  
  # RWD 
  f[['rwd']][,i] <- rwf(train.yi, drift = TRUE, h = hor)$mean
  fL[['rwd']][,i] <- rwf(train.yi, drift = TRUE, h = hor)$lower[,2]
  fH[['rwd']][,i] <- rwf(train.yi, drift = TRUE, h = hor)$upper[,2]
  
  
  # VAR
  lag.sel <- VARselect(var.data[1:TT, ])$selection[criterion]
  m.var <- VAR(var.data[1:TT, ], p= lag.sel)
  f[['var']][,i] <- predict(m.var, n.ahead = hor)$fcst[[1]][, 1]
  fL[['var']][,i] <- predict(m.var, n.ahead = hor)$fcst[[1]][, 2]
  fH[['var']][,i] <- predict(m.var, n.ahead = hor)$fcst[[1]][, 3]
  
  # LASSO VAR
  nlg = lag.sel
  mod.lasso <-constructModel(var.data[1:TT, ], p=lag.sel, 
                             "Basic", gran = c(500,10),
                             cv="Rolling",
                             MN=TRUE,
                             C = c(1, 1, 1,0, 1),
                             verbose=FALSE)
  cv.res = cv.BigVAR(mod.lasso)
  res =BigVAR.est(mod.lasso)
  tmp.Z = matrix(1, nn*nlg+1, 1);
  for (ik in 1:nlg) {
    tmp.Z[((ik-1)*nn+2):(ik*nn+1), ] <- var.data[TT-ik+1, ]
  }
  BB <- res$B[,,cv.res@index]
  for (kk in 1:hor) {
    f.t <- BB%*%tmp.Z
    f[['var_lasso']][kk, i] <- f.t[1,1]
    tmp.Z[(nn+2):((nlg)*nn+1), ] <- tmp.Z[2:((nlg-1)*nn+1), ]
    tmp.Z[2:(nn+1), ] <- f.t
  }
  
  #VECM
  ca.jo(var.data)
  f[['vecm']][, i]<-predict(vec2var(ca.jo(var.data[1:TT,]), r=1), n.ahead=hor)$fcst[[1]][,1]
  # fL[['vecm']][, i]<-predict(vec2var(ca.jo(var.data), r=1), n.ahead=hor)$lower[[1]][,1]
  # fH[['vecm']][, i]<-predict(vec2var(ca.jo(var.data), r=1), n.ahead=hor)$upper[[1]][,1]
  
  
  # Eval errors
  for (kk in 1:length(models)) {
    e[[kk]][, i] <- test.yi - f[[kk]][, i]
  }
  
  # Next sample
  TT <- TT+1
}

rmse <- matrix(NA, hor, length(models))
colnames(rmse) <- models

ratio <- matrix(NA, hor, length(models)-1)
colnames(ratio) <- models[1:(length(models)-1)]
for (kk in 1:length(models)) {
  tmp <- e[[kk]]*e[[kk]]
  rmse[, kk] <- sqrt(rowMeans(tmp))
}

for (kk in 1:length(models)-1) {
  ratio[, kk] <- rmse[, kk]/rmse[, ncol(rmse)]
}

rmse
ratio

#Checkresiduals
checkresiduals(m.arima)
checkresiduals(m.arimax)
checkresiduals(seasD)
checkresiduals(m.ets)

#Casuality тип а чтобы и нет то?
causality(m.var, cause = c('Weighted_rate', 'USDRUB', 'KeyRate', 'Inflation'))$Granger

#Hair Plots forecasting (cross validation)
for (m in models) {
  xx <- ts(var.data[, 1], start = c(2010, 1), frequency = 12)
  gg <- autoplot(xx, col = 'black')
  ddates <- time(xx)
  TT   <- nrow(train.var) + 1
  for (i in 1:ncol(f[[m]])) {
    xx0 <- ts(f[[m]][, i], start = ddates[TT], frequency = 12)
    gg <- gg + autolayer(xx0)
    TT <- TT+1
  }
  #print(gg)
  ggsave(filename = paste0('hairplot_forecast_', m, '.png'), plot = gg)
}

##Forecasting plots

for (m in models) {
  xx <- ts(var.data[, 1], start = c(2010, 1), frequency = 12)
  ddates <- time(xx)
  for (hh in 1:1) {
    TT   <- nrow(train.var) + hh
    xx0 <- ts(f[[m]][hh, ], start = ddates[TT], frequency = 12)
    xxl <- ts(fL[[m]][hh, ], start = ddates[TT], frequency = 12)
    xxh <- ts(fH[[m]][hh, ], start = ddates[TT], frequency = 12)
    gg <- autoplot(xx, col = 'black')
    gg <- gg + autolayer(xx0, color= 'red') +autolayer(xxl, color = 'blue', linetype= 'dashed')+
      autolayer(xxh, color = 'blue', linetype= 'dashed')
  #gg 
    ggsave(filename = paste0('h_fore_', m, '_h_', hh, '.png'),  
           plot = gg)
  }
}

