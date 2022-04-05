# #data <- ##SIMULATE NORMAL WITH CHANGE IN MEAN OR ELSE BROWNIAN MOTION WITH SOME SHIFT
# data <- c(rnorm(100, 0, 1), rnorm(80, 0, 3),rnorm(100, 0, 1), seq(1,4, length.out=120)+rnorm(120, 4, 1),
#           rnorm(100, 0, 1), rnorm(50, 8, 1), rnorm(500, 2, 1),rnorm(50, 6, 1), rnorm(150, 0, 1), rnorm(50, 6, 1))
# # # # 
# DisneyShares <- read.csv("D:/Zorie/Bristol/#Year3/MathsProject/Examples/DisneySharePrice.csv")
# Open_Disney <- DisneyShares$Open # The data records the share price at the opening and the closing of trade,
# data <- Open_Disney[1:7000]
# acf_object <- acf(Open_Disney,lag.max=14000) # This gives the autocorrelation function plot - given that lower lags are not corrected for at when calculating the autocorrelation at higher lags, most lags are seen as significant.
# pacf_object <- pacf(Open_Disney) # The partial autocorrelation function is our way around this: you should find that this gives that only the first and second lags are significant
# Disney_Today <- Open_Disney[3:length(Open_Disney)] # Disney's share price, not including the first two days.
# Disney_Yesterday <- Open_Disney[2:(length(Open_Disney)-1)] # Disney's share price, not including the first day, or the last day.
# Disney_TDBY <- Open_Disney[1:(length(Open_Disney)-2)] # Disney's share price, not including the last two days.
# Adjusted_Disney <-Disney_Today - pacf_object$acf[1]*Disney_Yesterday - pacf_object$acf[2]*Disney_TDBY # Adjusting the Disney share price according to values found from the partial autocorrelation function plot.
# data <- Adjusted_Disney[1:5000]
# data <- Open_Disney
# 
# # TeslaShares <- read.csv("D:/Zorie/Bristol/#Year3/MathsProject/Code/SharePrices/TSLA.csv")
# Open_Tesla <- TeslaShares$Open
data <- Open_Tesla[2200:3000]

plot(data, type="l")

#################
#CREATE G VECTOR#  
#################

##g is the pmf for distance between changepoints
##G_k is the sum up to that point, finding values up to N is overkill.
getG <- function(N, p){
  g <- dgeom(1:N, p)
  # plot(g, type="l")
  
  G <- rep(0,N)
  for (i in 1:N){
    G[i] <- sum(g[1:i])
  }
  # plot(G, type="l")
  # lines(pgeom(1:N, p))
  return(G)
}
# max(getG(300, 0.08))


##calculates probability of time between changepoints. 
calcP <- function(i){
  p <- rep(0,i) 
  if (i == 1){
    p[1] <- 1
  } else {
    for (k in 1:(i-1)){ # checking for Go, not defined = 0
      if (k==i-1){
        p[k] <- (G[i-k])
      } else {
        p[k] <- (G[i-k] - G[i-k-1])/(1 - G[i-k-1])
      }
    }
  }
  return(p)
}#need to pass i as i is used in other places. 

###Updating distance between changepoints in beta distribution.
getBG <- function(n, cpinterval){

  mean <- mean(cpinterval)   #### MEAN VS MEDIAN
  if (length(cpinterval) < 6){
    return(getG(N, 1/mean))
  }
  median <- median(cpinterval)
  target <- mean
  var <- var(cpinterval)
  if (is.na(var)){ #If not enough changepoints, it cannot yet find a variance
    print("var not valid")
    return(getG(N, 1/mean))
  }
  
  #Finding interval in which to scale
  mininterval <- min(cpinterval)
  maxinterval <- max(cpinterval)
  bmean <- (mean-mininterval)/(maxinterval-mininterval)
  bvar <- (var)/((maxinterval-mininterval)^2)
  
  #Checks Method Of Moments estimate is valid
  if (!(bvar<(bmean*(1-bmean)))){
    print("not good for mom estimate")
    return(getG(N, 1/mean))
  }
  
  a <- bmean * (((bmean*(1 - bmean))/bvar) - 1)
  b <- (1 - bmean)*(((bmean*(1 - bmean))/bvar) - 1)
  # Cumilative distribution function
  centre <- qbeta(0.5, a, b)
  deslen = N
  scaledlen = mean/centre
  upto <- 1
  outlen <- scaledlen
  if (scaledlen > deslen){
    upto <- pnorm(deslen/scaledlen, a, b)
    outlen <- deslen
  }
  by2 <- pbeta(seq(0, upto, length.out=outlen), shape1 = a, shape2 = b)
  if (outlen < deslen){
    by2 <- append(by2, rep(max(by2), (deslen-outlen)))
  }
  par(mfrow=c(1,1))
  print("returning BETA G")
  print(cpinterval)
  print(mean)
  print(qbeta(0.5, shape1 = a, shape2 = b)*scaledlen)
  plot(by2)
  return(by2)
}

getBetaG <- function(N, cpinterval) {
  mean <- mean(cpinterval)
  var <- var(cpinterval)
  if (is.na(var)){
    return(getG(N, 1/mean))
  }
  mininterval <- min(cpinterval) #minchange
  maxinterval <- max(cpinterval) #max(diff between changepoints)
  bmean <- (mean-mininterval)/(maxinterval-mininterval)
  bvar <- (var)/((maxinterval-mininterval)^2)
  
  if (!(bvar<(bmean*(1-bmean)))){
    return(getG(N, 1/mean))
  }
  a <- bmean * (((bmean*(1 - bmean))/bvar) - 1)
  b <- (1 - bmean)*(((bmean*(1 - bmean))/bvar) - 1)
  
  # Cumilative distribution function
  centre <- qbeta(0.5, a, b)
  print("BETA DISTRIBUTION!")
  print(cpinterval)
  
  centredmean <- mean/centre
  outlen <- min((mean/centre), N)
  if (centredmean > N){
    print("array will be too long")
    upto <- pbeta(centredmean/N, a, b)
    print(centredmean)
    print(upto)
  } else { upto <- 1 }
  print(mean)
  print(pbeta(N, shape1 = a, shape2 = b))
  G<- pbeta(seq(0, upto, length.out=N), shape1 = a, shape2 = b) 
  # if (outlen != N){
  #   G <- append(G, rep(1,(N-outlen)))
  # }
  #length chosen so there is a 50% chance of being above/below the mean
  plot(G)
  # print(G)
  return(G)
}
# getBetaG(100, 0.01^2, 0.01)
# getBetaG(200, 100^2, (100-1)/100^2)
# data[300:307]

##checks whether there has been a changepoint
isChangepoint <- function(row, varcount){
  len <- length(row)
  if (len >= minChange){ ##no changepoints before this
    # print("CHECKING FOR CHANGEPOINTS!")
    # print(len)
    # print(row)
    # print(n)
    # print(row[len])
    if ((row[len])>0.95){
      return(TRUE)
    }
    #return(varCheck(varcount))
  }
  return(FALSE)
}

###Finds w (probability of data being more extreme than current data point based on previous data point)
findw <- function(mean, sd){
  w <- pnorm(data[n], mean, sd)
  #if (data == NA){}
  if (w < 0.5){
    w <- 2*w
  } else {
    w <- 2*(1-w)
  }
  return(w)
}

calcMean <- function(i, curMean){
  newMean <- curMean*(i-1)
  newMean <- (newMean + data[n])/(i)
  return(newMean)
}

#using median absolute deviation to estimate the standard deviation
stdvMAD <- function(i, mean){
  if (i==1){  #just dont.
    return(curdev)
  }
  else{
    x <- abs(data[n]-mean)
    deviation <<- c(deviation[deviation < x], x, deviation[deviation >= x])
    mad <- median(deviation)
    return(1.4826*mad)
    # if (i < minChange){
    #   return(curdev)
    # }
    # else{
    #   return(1.4826*mad)
    # }
  }
}

##checking for change in variance by sorting into bins. 
############
varChange <- function(mean, dev){
  val = data[n]
  if (val > mean + dev){
    varcount[4] = varcount[4] + 1
  }
  if ((mean - dev <= val) & (val < mean)){
    varcount[2] = varcount[2] + 1
  }
  if ((mean < val) && (val < mean + dev)){
    varcount[3] = varcount[3] + 1
  }
  if (val < mean - dev){
    varcount[1] = varcount[1] + 1
  }
  return(varcount)
}

#checking if distribution of data is as expected.
##################
# varCheck <- function(varcount){
#   sum <- sum(varcount)
#   # print(varcount)
#   # print("in varCheck")
#   # print(timesincechangepoint)
#   # print(varcount)
#   # print(timesincechangepoint)
#   # print((varcount[1]+varcount[4])/sum)
#   # print(0.34 + (1/timesincechangepoint))
#   if (((varcount[1]+varcount[4])/sum) > (0.4 + (1/timesincechangepoint))){
#     print("variance has increased")
#     return(TRUE)
#   }
#   if (((varcount[2]+varcount[3])/sum) > (0.75 + (1/timesincechangepoint))) {
#     print("variance has decreased")
#     return(TRUE)
#   }
#   # if (((varcount[1]+varcount[2])/sum) > (0.6 + (1/timesincechangepoint))) {
#   #   print("trending down")
#   #   return(TRUE)
#   # }
#   # if (((varcount[3]+varcount[4])/sum) > (0.6 + (1/timesincechangepoint))) {
#   #   print("trending up")
#   #   return(TRUE)
#   # }
#   return(FALSE)
# }
#############

##takes current vector, and returns updated next one, which is longer by 1
updateModel <- function(row, mean=data[n], stdv){
  i <- length(row) # symbolises time since last changepoint. 
  newrow <- rep(0, i+1)
  for (j in 1:i){
    w <- findw(mean=rollingmean, sd=stdv)
    if (i == j){ ## as G(O) Not defined in R, = 0
      newrow[j] <- w * (1- G[i-(j-1)])*row[j]
    } else{
      newrow[j] =  w * (1- G[i-(j-1)]/(1-G[i-(j-1)-1]))* row[j]
    }
  }

  p <- calcP(i)
  newrow[i+1] <- p%*%row
  
  #NORMALISING
  newrow <- newrow/sum(newrow)
  return(newrow)
}


#################################
#INITIALISING RELEVANT VARIBALES#
#################################
# change_estimate <- 100
N <- length(data)
# G <- getBetaG(N, a, b)
G <- getG(N, 0.01)
minChange <- 20 #Min length between changepoints
changepoints <- rep(0,20) ##ASSUMED MAX NUMBER OF CHANGEPOINTS, AS APPENDING COPIES ENTIRE ROW ##Needs adjustments
numberchangepoints <- 1
recent_changepoint <- 1
row <- rep(1,1)
rollingmean <- 0
deviation <- c(NULL)
curdev <- 0.000001
timesincechangepoint <- 1
no_bins <- 4
varcount <- rep(0,no_bins)
obsGap <- 0.01
cpinterval <- NULL
rollingmeanline <- rep(0,N)
deviationline <- rep(0,N)
################################ 

###MAIN
for (n in 1:(N-1)){
  # print(n)
  i <- length(row)
  #print(i)
  rollingmean <- calcMean(i, rollingmean) #has to be here so variable is updated.
  rollingmeanline[n] <- rollingmean
  stdv <- stdvMAD(i, rollingmean)
  deviationline[n] <- stdv
  row <- updateModel(row, rollingmean, stdv)
  varcount <- varChange(rollingmean, stdv)

  timesincechangepoint <- timesincechangepoint +1
  if (isChangepoint(row, varcount)){
    print("FOUND CHANGEPOINT!!!!!")
    print(n)
    print(rollingmean)
    print(stdv)
    changepoints[numberchangepoints] <- n #####NEED TO CHANGE POSITION IN DATA
    numberchangepoints <- numberchangepoints + 1
    recent_changepoint <- n
    cpinterval <- append(cpinterval, timesincechangepoint) ##inefficient
    timesincechangepoint <- 1
    row <- 1
    rollingmean <- 0.000001
    curdev <- stdv
    deviation <- NULL
    # G <- getG(N, numberchangepoints/n)
    G <- getBG(N, cpinterval) #getBetaG(N, cpinterval)##
    varcount <- rep(0,no_bins)
  }
}


##########################
### OUTPUTS AND GRAPHS ###
##########################
print(changepoints)
plot(data, type="l")



# #DISNEY
# plot(data,xaxt="n",main="Disney",xlab="Time",ylab="Share Price",type="l") 
# labDates <- c("1970","1980","1990","2000","2010","2020")
# axis(side=1,at=c(1988,4514,7042,9570,12085,14601),labels=labDates)
# abline(v=changepoints, col="red")

# #TESLA
# plot(data,xaxt="n", main="Tesla",xlab="Time",ylab="Share Price",type="l") 
# labDates <- c("2012","2014","2016","2018","2020","2022")
# axis(side=1,at=c(384,886,1390,1893,2396,2901),labels=labDates)

# #TESLA
plot(data,xaxt="n", main="Tesla",xlab="Time",ylab="Share Price",type="l") 
labDates <- c("2020","2021","2022")
axis(side=1,at=c(196,449,701),labels=labDates)


lines(rollingmeanline, type = "l", col="blue")
lines((rollingmeanline + deviationline), type="l", col="green")
lines((rollingmeanline - deviationline), type="l", col="green")
# lines(data, type = "l")
abline(v=changepoints, col="red")
# locator()

print("done")

# deviationline
# plot(Open_Disney[5000:10000], type="l")
# abline(v=changepoints, col="red")
data[2330:2350]
############
  #AT THE END#
############

###### TO DO

######## Potential improvements. CHANGE IN MEAN!
### W BROWNIAN, 
#estimating sigma MAD (meadian absolute deviation, robust statistics) sequential hypothesis testing book, 
#####likelihood ratio for changing ratio
#####change in variance
# what data in bins, and see whether equal amounts fall in each bin. If in outer bins, variance grown, if in middle only

# optim
# geometric is a posterior estimator of gap between changepoints. 
# simulated annealing - choosing parameters, threshold parameter such that detect as many things as possible, with as few mistakes as possible. 
# - heuristic optimisation - looking for correct parameters and jump around, but as you get closer, you slowly converge to what you think is correct. 

###LOW PRIORITY
#starting again after x timepoint so matrix doesnt get ludicrously large? only if necessary
#consider splitting long times of no changepoints, keep 10 a couple so can detect immediate changepoints
#pointless efficiency improvement, lookup table for P, Too much space added to be worth it

#done
#get it to detect the changepoint peak
#make the matrix a vector
#MEAN FROM PREV values,
##Median absolute deviation for standard deviation
#### Learn the amount of time between changepoints rather than assuming 
# SAFEBAYES (could be bayesian about it)
# MCMC?
# pick a distribution, conjucgate prior of gemoetric, is beta, update each time a changepoint is found, 