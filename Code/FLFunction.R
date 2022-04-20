

calcP <- function(i, G){
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
}

getG <- function(N, p){
  g <- dgeom(1:N, p)

  G <- rep(0,N)
  for (i in 1:N){
    G[i] <- sum(g[1:i])
  }
  return(G)
}

getBetaG <- function(N, cpinterval){
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
  BetaG <- pbeta(seq(0, upto, length.out=outlen), shape1 = a, shape2 = b)
  if (outlen < deslen){
    BetaG <- append(BetaG, rep(max(BetaG), (deslen-outlen)))
  }
  par(mfrow=c(1,1))
  # print("returning BETA G")
  # print(cpinterval)
  # print(mean)
  # print(qbeta(0.5, shape1 = a, shape2 = b)*scaledlen)
  # plot(BetaG)
  return(BetaG)
}

isChangepoint <- function(row, threshold, minchange=1){
  len <- length(row)
  if (len >= minchange){ ##no changepoints before this
    # print("CHECKING FOR CHANGEPOINTS!")
    # print(len)
    # print(row)
    # print(n)
    # print(row[len])
    if ((row[len])>threshold){
      return(TRUE)
    }
  }
  return(FALSE)
}

findw <- function(datan, mean, sd){
  w <- pnorm(datan, mean, sd)
  if (w < 0.5){
    w <- 2*w
  } else {
    w <- 2*(1-w)
  }
  return(w)
}

updateModel <- function(row, datan, mean, stdv, G, i){
  #i <- length(row) # symbolises time since last changepoint. 
  newrow <- rep(0, i+1)
  for (j in 1:i){
    w <- findw(datan, mean, stdv)
    if (i == j){ ## as G(O) Not defined in R, = 0
      newrow[j] <- w * (1- G[i-(j-1)])*row[j]
    } else{
      newrow[j] =  w * (1- G[i-(j-1)]/(1-G[i-(j-1)-1]))* row[j]
    }
  }
  p <- calcP(i, G)
  newrow[i+1] <- p%*%row
  
  #NORMALISING
  newrow <- newrow/sum(newrow)
  return(newrow)
}

calcMean <- function(i, curMean, datan){
  newMean <- curMean*(i-1)
  newMean <- (newMean + datan)/(i)
  return(newMean)
}

FL <- function(data, threshold=0.97, stdv=1, est_interval=100, beta=FALSE, rollingmean=FALSE, MAD=FALSE, 
               minchange=20){
  #################################
  #INITIALISING RELEVANT VARIBALES#
  N <- length(data)
  changepoints <- 0
  numberchangepoints <- 1
  recent_changepoint <- -1
  row <- 1
  rollingmean <- 0
  deviation <- c(NULL)
  timesincechangepoint <- 1
  cpinterval <- NULL
  rollingmeanline <- rep(0,N)
  deviationline <- rep(0,N)
  ################################ 
  
  G <- getG(N, 1/est_interval)

  ###MAIN
  for (n in 1:(N)){
    print("start of loop")
    print(n)
    print(data[n])
    print(rollingmean)
    print(stdv)
    print(timesincechangepoint)
    datan <- data[n]
    row <- updateModel(row, datan, rollingmean, stdv, G, timesincechangepoint)
    
    if (isChangepoint(row, threshold, minchange)){
      print("CHANGEPOINT DETECTED")
      # print(n)
      changepoints[numberchangepoints] <- n #####NEED TO CHANGE POSITION IN DATA
      numberchangepoints <- numberchangepoints + 1
      recent_changepoint <- n
      deviation <- NULL
      # G <- getG(N, numberofchangepoints/n)
      if (beta){
        cpinterval <- append(cpinterval, timesincechangepoint) ##inefficient
        G <- getBetaG((N-n), cpinterval) #getBetaG(N, cpinterval)##
      }
      row <- 1
      rollingmean <- data[n]
      # curdev <- stdv
      timesincechangepoint <- 1
    } 
    else{ # UPDATE VARIABLES
      if (rollingmean){
        rollingmean <- calcMean(timesincechangepoint, rollingmean, data[n]) #has to be here so variable is updated.
        rollingmeanline[n] <- rollingmean
      } else {
        rollingmean <- data[n]
      }
      if (MAD){
        print(timesincechangepoint)
        if (timesincechangepoint==1){  #just dont.
        }
        else{
          x <- abs(data[n]-rollingmean)
          deviation <- c(deviation[deviation < x], x, deviation[deviation >= x])
          stdv <- median(deviation)*1.4826
          deviationline[n] <- stdv
        }
      }
      timesincechangepoint <- timesincechangepoint +1
    }
  }
  
  rollingmeanline <<- rollingmeanline
  deviationline <<- deviationline
  return(changepoints)
}

# input <- c(rnorm(100, 0, 1), rnorm(80, 0, 3),rnorm(100, 0, 1), seq(1,4, length.out=120)+rnorm(120, 4, 1),
#            rnorm(100, 0, 1), rnorm(50, 8, 1))
input <- c(rnorm(100, 0, 1), rnorm(80, 4, 1),rnorm(100, 0, 1), seq(1,4, length.out=120)+rnorm(120, 4, 1),
           rnorm(100, 0, 1), rnorm(50, 8, 1), rnorm(100, 0, 1), rnorm(80, 4, 1),rnorm(100, 0, 1))

# DisneyShares <- read.csv("D:/Zorie/Bristol/#Year3/MathsProject/Examples/DisneySharePrice.csv")
# Open_Disney <- DisneyShares$Open # The data records the share price at the opening and the closing of trade,
# input <- Open_Disney[1:400]

# TeslaShares <- read.csv("D:/Zorie/Bristol/#Year3/MathsProject/Code/SharePrices/TSLA.csv")
# Open_Tesla <- TeslaShares$Open
# input <- Open_Tesla[2200:3000]
# print(input[105:115])

changepointsOG <- FL(input)
changepointsBeta <- FL(input, beta=TRUE, rollingmean=TRUE, MAD=TRUE)


print(changepointsOG)
print(changepointsBeta)

plot(input, type="l")
lines(rollingmeanline, type = "l", col="blue")
lines((rollingmeanline + deviationline), type="l", col="green")
lines((rollingmeanline - deviationline), type="l", col="green")

abline(v=changepointsBeta, col="black")
abline(v=changepointsOG, col="blue")
changepointsjoint <- intersect(changepointsBeta, changepointsOG)
abline(v=changepointsjoint, col="red")
