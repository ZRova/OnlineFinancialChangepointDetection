
####TESTING RANDOM R features and syntax


#multiplying sections of lists
p <- 1:20
print(p)

q <- 1:5
print(q)

print(p*q)
print(p[2:6]*q)
print(p[1:4])

#updating variables in functions
x <- c(4, 6, 9, 10, 13, 15, 19)
x
i <- 11
x <- c(x[x < i], i, x[x >= i])
x
 

##getting a trending line with noise
data <- seq(1,5,)
print(data)
plot(data, type = "l")

returns <- function(x, y){
  return(x+y, x*y)
}

# a b <- returns(2,3)
# print(a)
## no must return list.

qnorm(0.2, 0, 1)



###Beta distribution in R
x <- seq(0, 1, length.out = 21)
print(x)
dbeta(x, 1, 1)
pbeta(x, 1, 1)
plot(dbeta(x,2,2))
plot(pbeta(x,2,2))

## Visualization, including limit cases:
pl.beta <- function(a,b, asp = if(isLim) 1, ylim = if(isLim) c(0,1.1)) {
  if(isLim <- a == 0 || b == 0 || a == Inf || b == Inf) {
    eps <- 1e-10
    x <- c(0, eps, (1:7)/16, 1/2+c(-eps,0,eps), (9:15)/16, 1-eps, 1)
  } else {
    x <- seq(0, 1, length.out = 1025)
  }
  fx <- cbind(dbeta(x, a,b), pbeta(x, a,b), qbeta(x, a,b))
  f <- fx; f[fx == Inf] <- 1e100
  matplot(x, f, ylab="", type="l", ylim=ylim, asp=asp,
          main = sprintf("[dpq]beta(x, a=%g, b=%g)", a,b))
  abline(0,1,     col="gray", lty=3)
  abline(h = 0:1, col="gray", lty=3)
  legend("top", paste0(c("d","p","q"), "beta(x, a,b)"),
         col=1:3, lty=1:3, bty = "n")
  invisible(cbind(x, fx))
}
pl.beta(3,1)

pl.beta(2, 4)
pl.beta(3, 7)
pl.beta(3, 7, asp=1)

pl.beta(0, 0)   ## point masses at  {0, 1}

pl.beta(0, 2)   ## point mass at 0 ; the same as
pl.beta(1, Inf)

pl.beta(Inf, 2) ## point mass at 1 ; the same as
pl.beta(3, 0)

pl.beta(Inf, Inf)# point mass at 1/2

# is += a thing in R
e <- 1
e
e = +1
e
#Conclusion. No. 

f <- 1
e&&f



varcount <- rep(0,4)
varcount[1] = 2
varcount[4] = 4
varcount
sum(varcount)
if ((varcount[1] + varcount[4])/10 > 0.6){
  print("greater")
} else {
  print("not")
}

min(4,5)

###Testing passing blank variable
blankFunction <- function(var1, var2 = TRUE){
  if (var1 & var2) {
    return("both true")
  }
  else return("not both true")
}
blankFunction(TRUE)
blankFunction(FALSE, FALSE)
blankFunction(FALSE)
blankFunction(TRUE, TRUE)


##Testing how well global variables work within a function. 
# 
# globalFunction <- function(data, var1, var2){
#   
# }


select <- seq(1:10)
select
select <- select[-1]
select



# Creating the Sequence
gfg = seq(0,1, by=0.1)

# Case 3
plot(gfg, dbeta(gfg, 2,2), xlab = "X",
     ylab = "Beta Density", type = "l",
     col = "Red")


# The Beta Distribution
plr.data <- data.frame(
  player_avg <- c(seq(0, 1, length.out=1000)),
  stringsAsFactors = FALSE
)

# Print the data frame.           
print(plr.data)
print(plr.data$player_avg)
by1 <- dbeta(plr.data$player_avg, shape1 = 5, shape2 = 8)
par(mar = rep(2,4))
plot(by1)


##################
###THIS IS HOW THE BETA DISTRIBUTION WORKS!!!###
######

mininterval <- 80 #minchange
maxinterval <- 150 #max(diff between changepoints)

mean <- 100 #seq(mininterval,maxinterval)
print(mean)
var <- 30

# array <- c(26, 40, 94, 205)
array <- c(26,  40, 94, 205, 130, 17,70,  63,  24,  45, 333,  62, 203, 186, 431,  18,  58,  46, 174,  26,  14,  62, 311,  28,  14, 183, 106,  35,  24, 219, 108, 289,
           409, 132, 173,  60, 193, 213,  38)
# array <- c(10,  15,  28,  15, 108,  10,  10,  32,  10)
mean <- mean(array)
median <- median(array)
var <- var(array)
mininterval <- min(array)
maxinterval <- max(array)

bmean <- (mean-mininterval)/(maxinterval-mininterval)
bvar <- (var)/((maxinterval-mininterval)^2)
a <- bmean * (((bmean*(1 - bmean))/bvar) - 1)
b <- (1 - bmean)*(((bmean*(1 - bmean))/bvar) - 1)
# Cummilative distribution function
centre <- qbeta(0.5, a, b)
centre
mean
maxlen = 600
mean/centre
upto <= 1
if (mean/centre > maxlen){
  upto <- pnorm((mean/centre)/maxlen, a, b)
}
outlen <- mean/centre
# if outlen > N 
pbeta(0.1, a, b)
qbeta(0.5, a, b)
upto
by2 <- pbeta(seq(0, upto, length.out=N), shape1 = a, shape2 = b)
par(mfrow=c(1,1))
max(by2)
plot(by2)
# locator()

par(mfrow=c(2,2))
plot(getG(1000, 0.01), main="First 1000 points of geom(0.01)")
plot(getG(100, 0.01), main="First 100 points of geom(0.01)")
plot(pbeta(seq(0, 1, length.out=1000), shape1 = a, shape2 = b), main="1000 points of beta")
plot(pbeta(seq(0, 1, length.out=100), shape1 = a, shape2 = b), main="1000 points of beta")

# Inverse Cummilative distribution function
by3 <- qbeta(seq(0, 1, length.out=maxinterval*4), shape1 = a, shape2 = b)
par(mar = rep(2,4))
plot(by3)

b4 <- rbeta(plr.data$player_avg, shape1 = a, shape2 = b)
par(mar = rep(2,4))
plot(density(b4), main = "Rbeta Plot")

array <- c(103, 180)
mean(array)
var(array)
#var((20))

naCheck <- NA
is.na(naCheck)
is.na(array)
!TRUE


###WHAT HAPPENS FOR VALUES OUTSIDE ARRAY
outside <- 1:50
outside
outside[25]
outside[60]

TRUE != FALSE

G <- pbeta(seq(0, 0.9, length.out = 100), shape1 = 3, shape2 = 2)
plot(G)
# locator()

pbeta(0.8, 2,2)
qbeta(0.8, 2,2)


# array <- c(26,  40, 94, 205, 130, 17,70,  63,  24,  45, 333,  62, 203, 186, 431,  18,  58,  46, 174,  26,  14,  62, 311,  28,  14, 183, 106,  35,  24, 219, 108, 289,
#            409, 132, 173,  60, 193, 213,  38)
# array <- c(218, 825,  20,  49, 200, 800, 50, 1500)
array <- c(105, 189, 110, 101,  51, 501,  22,  25, 156)
# array <- c(53,   32,   91,   48, 1301,  454, 20, 21)
mean <- mean(array)
median <- median(array)
var <- var(array)
mininterval <- min(array)
maxinterval <- max(array)

bmean <- (mean-mininterval)/(maxinterval-mininterval)
bvar <- (var)/((maxinterval-mininterval)^2)
a <- bmean * (((bmean*(1 - bmean))/bvar) - 1)
b <- (1 - bmean)*(((bmean*(1 - bmean))/bvar) - 1)
# Cummilative distribution function
centre <- qbeta(0.5, a, b)
# centre
deslen = 6000 # N
scaledlen = (mean)/centre
# scaledlen
# mean
# var
upto <- 1
outlen <- scaledlen
if (scaledlen > deslen){
  upto <- pbeta(deslen/scaledlen, a, b)
  outlen <- deslen
}
by2 <- pbeta(seq(0, upto, length.out=outlen), shape1 = a, shape2 = b)
if (outlen < deslen){
  by2 <- append(by2, rep(max(by2), (deslen-outlen)))
}
par(mfrow=c(1,1))
mean
scaledlen
qbeta(0.5, a, b)*scaledlen
# print(by2[mean])
# max(by2)
plot(by2)
# locator()

array <- 1:12
library(TTR)
array
rollno <- 2
rollup <- 2*rollno
EMA(array, rollno)
EMA(array, rollno)[rollno]
plot(array)
lines(EMA(array, rollno))

array[rollno:(2*rollno)]
array[3:6]

# mintimesincecp <- 10
# movingaverage <- function(curave, x, movavelen, timesincecp){
#   if moveavelen < mintimesincecp{
#     curave <- (curave*(timesincecp-1) + x)/timesincecp
#   }
#   else
#     cura
# }


####DEMO PLOTS

plot(pbeta(seq(0,1, length.out = 100), shape1 = a, shape2 = b))


