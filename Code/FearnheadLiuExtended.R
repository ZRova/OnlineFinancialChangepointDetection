#data <- ##SIMULATE NORMAL WITH CHANGE IN MEAN OR ELSE BROWNIAN MOTION WITH SOME SHIFT
data <- c(rnorm(80, 0, 1), rnorm(120, 8, 1), rnorm(100, 0, 1), rnorm(50, 8, 1))


plot(data, type="l")

#################
#CREATE G VECTOR#  
#################

##g is the pmf for distance between changepoints
##G_k is the sum up to that point
getG <- function(N){
  g <- dgeom(1:N, 0.01)
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

##checks whether there has been a changepoint
isChangepoint <- function(row){
  len <- length(row)
  if (len >= minChange){ ##no changepoints before this
    if ((row[len])>0.9){
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
G <- getG(N)
minChange <- 10 #Min length between changepoints
changepoints <- rep(0,10) ##ASSUMED MAX NUMBER OF CHANGEPOINTS, AS APPENDING COPIES ENTIRE ROW ##Needs adjustments
numberchangepoints <- 1
recent_changepoint <- 1
row <- rep(1,1)
rollingmean <- 0
deviation <- c(NULL)
curdev <- 1
print(deviation)
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
print("done")

############
#AT THE END#
############


###### TO DO
#starting again at each timepoint so matrix doesnt get ludicrously large

######## Potential improvements
### W BROWNIAN,  Median absolute deviation for standard deviation
#estimating sigma MAD (meadian absolute deviation, robust statistics) sequential hypothesis testing book, likelihood ratio for changing ratio
# what data in bins, and see whether equal amounts fall in each bin. If in outer bins, variance grown, if in middle only
#consider splitting long times of no changepoints, keep 10 a couple so can detect immediate changepoints
#remove the matrix, and replace with a vector

#pointless efficiency improvement, lookup table for P, Too much space added to be worth it

#### Learn the amount of time between changepoints rather than assuming 
# SAFEBAYES (could be bayesian about it)
# MCMC?
# pick a distribution, conjucgate prior of gemoetric, is beta, update each time a changepoint is found, 


#done
#get it to detect the changepoint peak
#make the matrix a vector
#MEAN FROM PREV values,