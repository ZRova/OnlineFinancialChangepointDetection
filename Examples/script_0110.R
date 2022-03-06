library('fdapace')
library('changepoint')

############################################
#SCRIPT FOR DEMONSTRATIONS - 1 OCTOBER 2021#
############################################

## Assessing Penalties/Model Fit via a Change-in-Mean example

data <- c(rnorm(100,0,1),rnorm(500,1,1),rnorm(400,-3,1)) # Creates data with two changes in mean under standard Gaussian noise

plot(1:1000,data,type="l") # Plots data

mean_all <- mean(data) # Computes MLE of mean across whole of the data series

abline(h=mean_all,col="red") # Places this mean on the plot

Squared_Residuals <- (data-mean_all)^2 # Creates a series of the squared residuals (i.e. squared distances of each point in original series from mean)

Null_RSS <- sum(Squared_Residuals) # Adds up the series of squared residuals created above in the no changepoint setting.
Null_RSS

plot(1:1000,data,type="l") # Re-plots the data

change_object <- cpt.mean(data,method="PELT",penalty="BIC") # Performs change-in-mean detection on the data using PELT with the Bayesian Information Criterion penalty.

abline(v=change_object@cpts,col="blue") # We fit the changepoints found by PELT.

##We adjust the series according to these discovered changes, in order to assess the model fit.

altered_series <- rep(0,1000)

altered_series[1:change_object@cpts[1]] <- data[1:change_object@cpts[1]]-change_object@param.est$mean[1]

plot(1:1000,altered_series,type="l")

altered_series[(change_object@cpts[1]+1):change_object@cpts[2]] <- data[(change_object@cpts[1]+1):change_object@cpts[2]] -change_object@param.est$mean[2] 

plot(1:1000,altered_series,type="l")

altered_series[(change_object@cpts[2]+1):change_object@cpts[3]] <- data[(change_object@cpts[2]+1):change_object@cpts[3]] -change_object@param.est$mean[3] 

plot(1:1000,altered_series,type="l")

RSS <- sum(altered_series^2) # We compute the residual sum of squares of the model implied by PELT (with the two changepoints). The result (should be!) much lower than the residual sum of squares found in the no changepoint setting.
RSS

Total_Cost <- RSS + 2*change_object@pen.value # We compute the total cost, which is the residual sum of squares plus the number of changepoints discovered multiplied by the penalty value.
Total_Cost

#####################

#Change in mean - what should the penalty be?

plot(1:1000,data,type="l") # Same example as above

change_object <- cpt.mean(data,method="PELT",penalty="CROPS",pen.value=c(1,10000)) # Running Changepoints for a Range Of PenaltieS (CROPS) - we look for minimum penalty which gives us a certain number of changes for penalty values between 1 and 1000.

change_mat <- change_object@pen.value.full # Outputs the minimum penalties in increasing order

plot(1:(length(change_mat)-2),change_mat[(length(change_mat)-1):2],type="l",xlab="Number of Changepoints",ylab="Minimum Penalty Required") # These two lines show the elbow plots we discussed - note that this one is dominated by the values when 'Number of Changepoints' is small, hence the re-plotting below.

plot(2:(length(change_mat)-2),change_mat[(length(change_mat)-2):2],type="l",xlab="Number of Changepoints",ylab="Minimum Penalty Required") # These two lines show the elbow plots we discussed.

# For this example, a penalty of anything above around 6 or 7 looks sensible from the elbow plot's perspective. Note that the BIC = 13.8 in our example.

###########################

DisneyShares <- read.csv("D:/Zorie/Bristol/#Year3/Project/DisneySharePrice.csv") # Disney Share Price - I've put this time series in the files section of the Teams channel as well.

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

###################
#ROC Curve Example#
###################

lots_of_changes <- c(rnorm(100,0,1),rnorm(100,1,1),rnorm(100,3,1),
                     rnorm(100,0,1),rnorm(100,-2,1),rnorm(100,4,1),
                     rnorm(100,2,1),rnorm(100,0,1),rnorm(100,6,1),
                     rnorm(100,3,1),rnorm(100,-1,1),rnorm(100,8,1),
                     rnorm(100,0,1),rnorm(100,1,1),rnorm(100,3,1),
                     rnorm(100,0,1),rnorm(100,-2,1),rnorm(100,4,1),
                     rnorm(100,2,1),rnorm(100,0,1),rnorm(100,6,1),
                     rnorm(100,3,1),rnorm(100,-1,1),rnorm(100,8,1),
                     rnorm(100,5,1),rnorm(100,0,1)) # A synthetic data example with 25 changepoints: there is a change in mean at every 100 points.

plot(1:2600,lots_of_changes,type="l") # Plot of the data above.

change_object<-cpt.mean(lots_of_changes,penalty="CROPS",method="PELT",pen.value=c(1,10000)) # CROPS plot to derive a series of "critical penalties" - we'll use these in the ROC curve (although we could equally just use any set of numbers we choose).

good_set <- c(seq(97,2497,100),seq(98,2498,100),seq(99,2499,100),seq(100,2500,100),seq(101,2501,100),seq(102,2502,100),seq(103,2503,100)) # If a change detected by PELT (under a given penalty) falls inside this "good set" (i.e. gives an error of at most three in finding a change).

howmanynotin <- function(list_of_changes,good_set){ # This function counts how many false alarms are raised.
  logical_vec <- rep(0,length(list_of_changes))
  for(i in 1:length(list_of_changes)){
    if(list_of_changes[i]%in%good_set){
      logical_vec[i] <- 1
    }
  }
  return(length(list_of_changes)-sum(logical_vec))
}

howmanyfound <- function(list_of_changes,true_locations){ # This function counts how many true changes are detected.
  num_real <- length(true_locations)
  found <- rep(0,num_real)
  for(i in 1:num_real){
    error_set <- seq(true_locations[i]-3,true_locations[i]+3,1)
    if(sum(error_set%in%list_of_changes)>0){
      found[i]<-1
    }
  }
  return(sum(found))
}

number_of_points <- dim(change_object@cpts.full)[1] # The number of penalties we are going to place in the ROC plot.

false_positive_rate <- rep(0,number_of_points) # Initialising the vector of false positives.
true_positive_rate <- rep(0,number_of_points) # Initialising the vector of true positives.

for(i in 1:number_of_points){ # This loop compute the false and true positive rates for each penalty.
  number<-howmanynotin(change_object@cpts.full[i,][!is.na(change_object@cpts.full[i,])],good_set)
  
  false_positive_rate[i] <- number/length(change_object@cpts.full[i,][!is.na(change_object@cpts.full[i,])])
  
  found<-howmanyfound(change_object@cpts.full[i,],seq(100,2500,100))
  
  true_positive_rate[i] <- found/25
}

plot(false_positive_rate,true_positive_rate,type="p",xlab="False Positive Rate",ylab="True Positive Rate",main="Example ROC Curve",pch=20) # This computes the ROC plot.
abline(v=0.05,col="red") # This line indicates the maximum value of the false positive rate we are prepared to tolerate (this may differ depending on the practitioner). 
abline(h=0.90,col="red") # This line indicates the minimum value of the true positive rate we are prepared to tolerate (this may differ depending on the practitioner).


#Change in variance example

mydata <- c(rnorm(100,0,1),rnorm(100,0,3))

plot(1:200,mydata,type="l")


results <- cpt.var(mydata,penalty="BIC",know.mean=TRUE,mu=0,method="AMOC",test.stat="Normal",class="TRUE")

plot(1:200,mydata,type="l")

abline(v=results@cpts,col="blue")

######################

########################
#Geometric Distribution#
########################

sequence <- 0:30

success_prob <- 0.20

geometric_density <- dgeom(sequence,success_prob)

plot(sequence+1,geometric_density,type="p",pch=20)

#################
#Brownian Motion#
#################

n <- 10000

W <- Wiener(n = 1, pts = seq(0, 1, length = n), sparsify = NULL, K = 500)

plot(1:length(W),W,type="l")
