#data <- ##SIMULATE NORMAL WITH CHANGE IN MEAN OR ELSE BROWNIAN MOTION WITH SOME SHIFT
#data <- c(rnorm(80, 0, 1), seq(1,4, length.out=120)+rnorm(120, 5, 1), rnorm(100, 0, 1), rnorm(50, 8, 1))

DisneyShares <- read.csv("D:/Zorie/Bristol/#Year3/MathsProject/Examples/DisneySharePrice.csv") # Disney Share Price - I've put this time series in the files section of the Teams channel as well.

Open_Disney <- DisneyShares$Open # The data records the share price at the opening and the closing of trade, as well as the high and low achieved by the share price each day. We focus on the opening price here. (Perhaps you could investigate whether a similar pattern holds for the closing, high and low prices?)
plot(Open_Disney[1:4000], type = "l")

data <- Open_Disney

# DisneyShares <- read.csv("D:/Zorie/Bristol/#Year3/MathsProject/Examples/DisneySharePrice.csv") # Disney Share Price - I've put this time series in the files section of the Teams channel as well.
# 
# data <- DisneyShares$Open
plot(data, type="l")

#################
#CREATE G VECTOR#  
#################

##g is the pmf for distance between changepoints
##G_k is the sum up to that point
getG <- function(N, p){
  g <- dgeom(1:N, p)
  #plot(g, type="l")
  
  G <- rep(0,N)
  for (i in 1:N){
    G[i] <- sum(g[1:i])
  }
  return(G)
}

##calculates probability of time between changepoints. 
calcP <- function(i){
  p <- rep(0,i) 
  if (i == 1){
    p[1] <- 1
  } else {
    for (k in 1:(i-1)){ # checking for Go, not defined = 0
      if (k==i-1){
        p[k] <- (G[i-k])
      }else{
        p[k] <- (G[i-k] - G[i-k-1])/(1 - G[i-k-1])
      }
    }
  }
  return(p)
}#need to pass i as i is used in other places. 

###Updating distance between changepoints in beta distribution.
# updateG <- function(changepoints, position){
#   betaNew <- 
# }

##checks whether there has been a changepoint
isChangepoint <- function(row){
  len <- length(row)
  if (len >= minChange){ ##no changepoints before this
    if ((row[len])>0.95){
      return(TRUE)
    }
  }
  return(FALSE)
}

###Finds w (probability of data being more extreme than current data point based on previous data point)
findw <- function(data, mean, sd){
  w <- pnorm(data, mean, sd)
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
    if (i < minChange){
      return(curdev)
    }
    else{
      return(1.4826*mad) 
    }
  }
}

##checking for change in variance by sorting into bins. 
varcount <- rep(0,10)
probRange <- qnorm(seq(0,1, by=0.1), mean, stdv)
changeVar <- rep(0,10)

##takes current vector, and returns updated next one, which is longer by 1
updateModel <- function(row, mean=data[n], stdv=1){
  i <- length(row) # symbolises time since last changepoint. 
  newrow <- rep(0, i+1)
  for (j in 1:i){
    w <- findw(data[n+1], mean=rollingmean, sd=stdv)
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
N <- length(data)
G <- getG(N, 0.01)
minChange <- 10 #Min length between changepoints
changepoints <- rep(0,20) ##ASSUMED MAX NUMBER OF CHANGEPOINTS, AS APPENDING COPIES ENTIRE ROW ##Needs adjustments
numberchangepoints <- 1
recent_changepoint <- 1
row <- rep(1,1)
rollingmean <- 0
deviation <- c(NULL)
curdev <- 1
################################ 

###MAIN
for (n in 1:(N-1)){
  i <- length(row)
  #print(i)
  rollingmean <- calcMean(i, rollingmean) #has to be here so variable is updated.
  stdv <- stdvMAD(i, rollingmean)
  row <- updateModel(row, rollingmean, stdv)
  #print(row)
  if (isChangepoint(row)){
    print("FOUND CHANGEPOINT!!!!!")
    print(n)
    print(stdv)
    print(rollingmean)
    changepoints[numberchangepoints] <- n #####NEED TO CHANGE POSITION IN DATA
    numberchangepoints <- numberchangepoints + 1
    recent_changepoint <- n
    row <- 1
    rollingmean <- 0
    deviation <- NULL
  }
}

print(changepoints)
plot(data[1:7000], type = "l")
abline(v=changepoints, col="red")
print("done")

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


#### Learn the amount of time between changepoints rather than assuming 
# SAFEBAYES (could be bayesian about it)
# MCMC?
# pick a distribution, conjucgate prior of gemoetric, is beta, update each time a changepoint is found, 


#done
#get it to detect the changepoint peak
#make the matrix a vector
#MEAN FROM PREV values,
##Median absolute deviation for standard deviation