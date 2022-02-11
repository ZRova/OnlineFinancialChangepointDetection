
#data <- ##SIMULATE NORMAL WITH CHANGE IN MEAN OR ELSE BROWNIAN MOTION WITH SOME SHIFT
data <- c(rnorm(8, 0, 1), rnorm(12, 8, 1))
N <- length(data)

plot(data, type="l")

#################
#CREATE G VECTOR#  
#################

g <- dgeom(1:N, 0.1)
#plot(g, type="l")

G <- rep(0,N)
for (i in 1:N){
  G[i] <- sum(g[1:i])
}
#plot(G, type="l")


##checks whether there has been a changepoint
isChangepoint <- function(row){
  len <- length(row)
  if (len >= 3){
    if ((row[len] + row[len-1] + row[len-2])>1.1){
      print(len)
      return(TRUE)
    }
  }
  return(FALSE)
}

###Finds w (probability of data being more extreme than current data point based on previous data point)
findw <- function(data, mean, sd){
  w <- pnorm(data, mean, sd)
  if (w < 0.5){
    w <- 2*w
  } else {
    w <- 2*(1-w)
  }
  print(w)
  return(w)
}

##takes current vector, and returns updated next one, which is longer by 1
updateModel <- function(row){
  i <- length(row)
  newrow <- rep(0, i+1)
  for (j in 1:i){
    print("data")
    print(data[i+1])
    w <- findw(data[i+1], data[i], sd=1)
    if (i == j){ ## as G(O) Not defined in R, = 0
      #output_matrix[i+1, j] <- w * (1- G[i-(j-1)])*output_matrix[i,j]  #################HEREEEEEE
      newrow[j] <- w * (1- G[i-(j-1)])*row[j]
    } else{
      #output_matrix[i+1,j] = w * (1- G[i-(j-1)]/(1-G[i-(j-1)-1]))*output_matrix[i,j]
      newrow[j] =  w * (1- G[i-(j-1)]/(1-G[i-(j-1)-1]))* row[j]
    }
  }
  
  ## J = i+1
  j <- i+1
  p <- rep(0,i) ##WHAT IS K
  if (i == 1){
    p[1] <- 1
  } else {
    for (k in 1:(i-1)){ # checking for Go, not defined = 0
      if (k==i-1){
        p[k] <- (G[i-k])
      }else{
        #print(G[i-k])
        p[k] <- (G[i-k] - G[i-k-1])/(1 - G[i-k-1])
      }
    }
  }
  
  #output_matrix[i+1,j] <- p%*%output_matrix[i,1:i] 
  newrow[j] <- p%*%row
  
  #NORMALISING
  #output_matrix[i+1,] <- output_matrix[i+1,]/sum(output_matrix[i+1,])
  newrow <- newrow/sum(newrow)
  return(newrow)
}

##ASSUMED MAX NUMBER OF CHANGEPOINTS, AS APPENDING COPIES ENTIRE ROW ##Needs adjustments
changepoints <- rep(0,10) 
numberchangepoints <- 0
recent_changepoint <- 0

##Initialising row
row <- rep(1,1)
print(row)

for (x in 1:(N-1)){
  ##START
  row <- updateModel(row)
  if (isChangepoint(row)){
    print("FOUND CHANGEPOINT!!!!!" + string(x))
    changepoints[numberchangepoints] = x #####NEED TO CHANGE POSITION IN DATA
    numberchangepoints = numberchangepoints + 1
    newdata <- data[i:(N-1)]
    data <- newdata
    row <- 1
  }
}

print(changepoints)
plot(row, type = "l")
print("done")

############
#AT THE END#
############

#plot(output_matrix[52,], type="l")

#After iterating through the above for each time point, you can then take a look
#at the posterior distribution for each time point. (For example, you could do a
#plot of some of the rows and see what you see, and perhaps turn this into an 
#animation). For time points after the change, do you see a spike around the 
#true changepoint in the posterior?

print("done")

###### TO DO
#get it to detect the changepoint peak
#make the matrix a vector
#starting again at each timepoint so matrix doesnt get ludicrously large

######## Potential improvements
### W BROWNIAN, MEAN FROM PREV values, Median absolute deviation for standard deviation
#estimating sigma MAD (meadian absolute deviation, robust statistics) sequential hypothesis testing book, likelihood ratio for changing ratio
# what data in bins, and see whether equal amounts fall in each bin. If in outer bins, variance grown, if in middle only
#consider splitting long times of no changepoints, keep 10 a couple so can detect immediate changepoints
#remove the matrix, and replace with a vector

#### Learn the amount of time between changepoints rather than assuming 
# SAFEBAYES (could be bayesian about it)
# MCMC?
# pick a distribution, conjucgate prior of gemoetric, is beta, update each time a changepoint is found, 
