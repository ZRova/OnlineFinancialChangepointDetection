#data <- ##SIMULATE NORMAL WITH CHANGE IN MEAN OR ELSE BROWNIAN MOTION WITH SOME SHIFT
data <- c(rnorm(8, 0, 1), rnorm(12, 8, 1), rnorm(10, 0, 1), rnorm(5, 8, 1))
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
  if (len >= 3){ ##no changepoints before this
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

##takes current vector, and returns updated next one, which is longer by 1
updateModel <- function(row){
  i <- length(row)
  newrow <- rep(0, i+1)
  for (j in 1:i){
    w <- findw(data[n+1], data[n], sd=1)
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
numberchangepoints <- 1
recent_changepoint <- 1

##Initialising row
row <- rep(1,1)
print(row)

for (n in 1:(N-1)){
  ##START
  row <- updateModel(row)
  print(row)
  if (isChangepoint(row)){
    print("FOUND CHANGEPOINT!!!!!")
    print(n)
    changepoints[numberchangepoints] = n #####NEED TO CHANGE POSITION IN DATA
    numberchangepoints = numberchangepoints + 1
    #print(data)
    #newdata <- rep(0,(N-x))
    #newdata <- data[x+1:(N)] #### CHANGE THIS TOO
    #data <- newdata
    #print("new data")
    #print(newdata)
    row <- 1
  }
}

print(changepoints)
#plot(row, type = "l")
print("done")

