library('fdapace')
library('changepoint')

DisneyShares <- read.csv("D:/Zorie/Bristol/#Year3/MathsProject/Examples/DisneySharePrice.csv") # Disney Share Price - I've put this time series in the files section of the Teams channel as well.

Open_Disney <- DisneyShares$Open # The data records the share price at the opening and the closing of trade, as well as the high and low achieved by the share price each day. We focus on the opening price here. (Perhaps you could investigate whether a similar pattern holds for the closing, high and low prices?)

plot(1:length(Open_Disney),Open_Disney,type="l",xlab="Time (Trading Days since January 1962)",ylab="Share Price ($)",main="Disney Share Price")

plot(1:length(Open_Disney),Open_Disney,xaxt="n",xlab="Time",ylab="Share Price",type="l") # These three lines re-plot the data to give a nice time-labelled x-axis.
labDates <- c("1970","1980","1990","2000","2010","2020")
axis(side=1,at=c(1988,4514,7042,9570,12085,14601),labels=labDates)

change_object<-cpt.mean(Open_Disney,penalty="BIC",method="PELT") # Performing changepoint detection on the opening share price data. (Change-in-mean using PELT with a BIC penalty -> you should find that this overfits and gives too many changepoints.)

abline(v=change_object@cpts,col="blue") # This plots the changepoints we have found.

acf_object <- acf(Open_Disney,lag.max=14000) # This gives the autocorrelation function plot - given that lower lags are not corrected for at when calculating the autocorrelation at higher lags, most lags are seen as significant.

pacf_object <- pacf(Open_Disney) # The partial autocorrelation function is our way around this: you should find that this gives that only the first and second lags are significant

Disney_Today <- Open_Disney[3:length(Open_Disney)] # Disney's share price, not including the first two days.
Disney_Yesterday <- Open_Disney[2:(length(Open_Disney)-1)] # Disney's share price, not including the first day, or the last day.
Disney_TDBY <- Open_Disney[1:(length(Open_Disney)-2)] # Disney's share price, not including the last two days.

Adjusted_Disney <-Disney_Today - pacf_object$acf[1]*Disney_Yesterday - pacf_object$acf[2]*Disney_TDBY # Adjusting the Disney share price according to values found from the partial autocorrelation function plot.

plot(1:length(Adjusted_Disney),Adjusted_Disney,type="l",xaxt="n",xlab="Time (Trading Days since January 1962)",ylab="Autocorrelation Corrected Share Price",main="Disney Share Price") # This series is a little more stationary than the raw data we plotted above.
labDates <- c("1970","1980","1990","2000","2010","2020")
axis(side=1,at=c(1988,4514,7042,9570,12085,14601),labels=labDates)

change_object <- cpt.meanvar(Adjusted_Disney,penalty="BIC",method="PELT") # Detecting changes in mean and variance.
abline(v=change_object@cpts,col="blue") # Again, we find a number of (i.e. too many) changes!

change_object <- cpt.meanvar(Adjusted_Disney,penalty="CROPS",method="PELT",pen.value=c(28.8,1000)) # Starting with the BIC as our "lower bound" for the penalty.
change_mat <- change_object@pen.value.full
no_changes <- change_object@cpts.full

##We now use the result of applying CROPS to find the "correct" penalty value.

xseq <- rep(0,length(change_mat)-2)
for(i in 1:(length(change_mat)-2)){
  xseq[i] <- length(no_changes[length(change_mat)-i,])-sum(is.na(change_object@cpts.full[length(change_mat)-i,]))+1
}

plot(xseq,change_mat[(length(change_mat)-1):2],type="l",xlab="Number of Changepoints",ylab="Minimum Penalty Required") # This gives us the elbow plot.

change_object <- cpt.meanvar(Adjusted_Disney,penalty="Manual",method="PELT",pen.value=200) # We use a penalty value of 200, as this corresponds to a sensible value from the elbow plot.

plot(1:length(Adjusted_Disney),Adjusted_Disney,type="l",xaxt="n",xlab="Time (Trading Days since January 1962)",ylab="Autocorrelation Corrected Share Price",main="Disney Share Price") # Replot the data to see the new changes.
labDates <- c("1970","1980","1990","2000","2010","2020")
axis(side=1,at=c(1988,4514,7042,9570,12085,14601),labels=labDates)

abline(v=change_object@cpts,col="blue") # Changepoints found using the manual penalty found from the elbow plot.

Log_Adjusted_Disney <- log(Disney_Today/Disney_Yesterday) # Performing the standard log-adjustment transform

plot(1:length(Log_Adjusted_Disney),Log_Adjusted_Disney,type="l",xaxt="n",xlab="Time (Trading Days since January 1962)",ylab="Autocorrelation Corrected Share Price",main="Disney Share Price")
labDates <- c("1970","1980","1990","2000","2010","2020")
axis(side=1,at=c(1988,4514,7042,9570,12085,14601),labels=labDates)

change_object <- cpt.var(Log_Adjusted_Disney,penalty="BIC",method="PELT") # Fitting changepoints to the log-transformed shares

abline(v=change_object@cpts,col="blue") # The model overfits again.

change_object <- cpt.var(Log_Adjusted_Disney,penalty="CROPS",method="PELT",pen.value=c(28.8,1000)) # Running CROPS again, with the minimum penalty again 28.8.
change_mat <- change_object@pen.value.full
no_changes <- change_object@cpts.full

xseq <- rep(0,length(change_mat)-2)
for(i in 1:(length(change_mat)-2)){
  xseq[i] <- length(no_changes[length(change_mat)-i,])-sum(is.na(change_object@cpts.full[length(change_mat)-i,]))+1
}

plot(xseq,change_mat[(length(change_mat)-1):2],type="l",xlab="Number of Changepoints",ylab="Minimum Penalty Required") # Elbow plot.

change_object <- cpt.var(Log_Adjusted_Disney,penalty="Manual",method="PELT",pen.value=150) # PELT with a manual penalty derived from the elbow plot.

plot(1:length(Log_Adjusted_Disney),Log_Adjusted_Disney,type="l",xaxt="n",xlab="Time (Trading Days since January 1962)",ylab="Autocorrelation Corrected Share Price",main="Disney Share Price")
labDates <- c("1970","1980","1990","2000","2010","2020")
axis(side=1,at=c(1988,4514,7042,9570,12085,14601),labels=labDates)

abline(v=change_object@cpts,col="blue") # Plotting the changes found above.

yw_object <- ar.yw(Open_Disney,order.max=10) # An alternative way of performing the stationarity correction involves computing the AR coefficients via Yule-Walker.

