################################################################################
#Fearnhead and Liu pseudo-code to compute the posterior of the most recent     #
#change location at each time point. This pseudo-code will detail the          #
#production of a matrix in which the ith row denotes P(C_i = .|y_{1:i}) for i  #
#from 1 to N, where N is the length of the inputted data.                      #
################################################################################

############
#DATA INPUT#
############

#Note we assume that we have seen all the data in advance for the purposes of 
#the R code. Therefore, we detect changes in an "as-if" online fashion.

#data <- ##SIMULATE NORMAL WITH CHANGE IN MEAN OR ELSE BROWNIAN MOTION WITH SOME SHIFT
data <- c(rnorm(50, 0, 1), rnorm(50, 8, 1))
#N <- length(data)
N <- length(data)

plot(data, type="l")

#################
#CREATE G VECTOR#
#################

#Create G vector in which the entries are G[0], ..., G[N], where G[k] is the G_k 
#defined in the general description. (I.e. use rgeom and sum.)

##g is the pmf for distance between changepoints
##G_k is the sum up to that point

g <- dgeom(1:N, 0.01)
plot(g, type="l")

G <- rep(0,N)
for (i in 1:N){
  G[i] <- sum(g[1:i])
}
plot(G, type="l")

################################
#INITIALISE FINAL OUTPUT MATRIX#
################################

#output_matrix <- matrix(rep(0,N*N),N,N) 
output_matrix <- matrix(rep(0,N*N),N,N) 
##This is a n*n matrix of 0s

#Note 1: recall that the ith row will be the posterior distribution for the most 
#recent changepoint after observing the ith timepoint. 

#Note 2: For ease of computation, output_matrix[i,j] = P(C_i = j-1|y_{1:i}), so
#for example output_matrix[1,1] = P(C_1 = 0|y_1) = 1, as C_1 can only take the
#value 0.

#Note 3: By Note 2 and the definition of the Markov Chain C_t, 
#output_matrix[i,j] = 0 if j > i.

###########
#FIRST ROW#
###########

#output_matrix[1,1] <- 1 #by Note 2 above.
output_matrix[1,1] <- 1

##Checking and other functions
isChangepoint <- function(row){
  #length of nonzero row
  len <- length(row)
  if (len>= 3){
    if ((row[len] + row[len-1] + row[len-2])>0.9){
      print(len)
      return(TRUE)
    }
  }
}

#################
#SUBSEQUENT ROWS#
#################
changePoints <- c(0)

#Suppose we have generated the rows of the matrix up to and including row i (for
#example by using a for loop). We now proceed inductively to generate row i+1.

##GENERATE output_matrix[i+1,j] for j = 1, ..., i.

#output_matrix[i+1,j] = w_{j-1,i+1} * (1 - G(i-(j-1))/(1-G(i-(j-1)-1))*output_matrix[i,j]

#Note that this w-value above is the p-value for how extreme the most recent 
#piece of data is given the last piece of data.
###### W = P(Y_T+1 | C_T+1 = j, Y_1:t) 
###### Y IS THE DATA!!!!!

##GENERATE output_matrix[i+1,j] for j = i+1

#First generate vector p = (G[i-k] - G[i-k-1])/(1 - G[i-k-1]), for k = 0, ..., i-1.

#output_matrix[i+1,j] = p*output_matrix[i,1:i]

##THEN normalise:

#output_matrix[i+1,] = output_matrix[i+1,]/sum(output_matrix[i+1,])

for (i in 1:(N-1)){
  for (j in 1:i){
    #w = probablity(data[i+1]| output_matrix[i, j], data[1:i]) (outer tails)
    w <- pnorm(data[i+1], mean=data[i], sd=1)
    if (w < 0.5){
      w <- w*2
    } else {
      w <- 2*(1-w)
    }
    #print(w)
    if (i == j){ ## as G(O) Not defined in R, = 0
      output_matrix[i+1, j] <- w * (1- G[i-(j-1)])*output_matrix[i,j]
    } else{
      output_matrix[i+1,j] = w * (1- G[i-(j-1)]/(1-G[i-(j-1)-1]))*output_matrix[i,j]
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
  
  output_matrix[i+1,j] <- p%*%output_matrix[i,1:i]
  
  #NORMALISING
  output_matrix[i+1,] <- output_matrix[i+1,]/sum(output_matrix[i+1,])
  
  ##### MATRIX COMPLETE
  
#  if (isChangepoint(output_matrix[i+1,])){
#    changePoints <- append(changePoints[i])
#    #AND RESTART
#  }
  
}


#print(output_matrix)
############
#AT THE END#
############

plot(output_matrix[52,], type="l")

#After iterating through the above for each time point, you can then take a look
#at the posterior distribution for each time point. (For example, you could do a
#plot of some of the rows and see what you see, and perhaps turn this into an 
#animation). For time points after the change, do you see a spike around the 
#true changepoint in the posterior?

library(animation)
library(ggplot2)
library(gganimate)

df <- as.data.frame(output_matrix)
colnames(df)<-c(1:N)

myplot <- ggplot(df, aes(x=x, y=y)) + geom_point()

animateplot <- myplot + transition_time(colnames(df)) + xlab("time") + ylab("prob of changepoint")

animate(animateplot)

print("leaving data")

# saveGIF({
# }
# )

desc = c("This is a super cool example of Gradient Descent")
saveGIF({
  f1 = function(x, y) x^2 + 3 * sin(y)
  xx = grad.desc(f1, pi * c(-2, -2, 2, 2), c(-2 * pi, 2))
  xx$persp(col = "lightblue", theta = 30, phi = 30)
},title = "Demo of Gradient Descent", description = desc, verbose = FALSE)

saveHTML({
  
})

print("done")

###### TO DO
#get it to detect the changepoint peak
#starting again at each timepoint so matrix doesnt get ludicrously large




######## Potential improvements
### W BROWNIAN, MEAN FROM PREV values, Median absolute deviation for standard deviation
#estimating sigma
#consider splitting long times of no changepoints, keep 10 a couple so can detect immediate changepoints
#remove the matrix, and replace with a vector






