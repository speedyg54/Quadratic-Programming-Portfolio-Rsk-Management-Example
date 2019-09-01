## Example in Quadratic Programming
# Assume we want to minimize: 1/2 x^*TC*x - (0 5 0)^T*x
# under the constraints:      A x <= b
# with b = (8,-2, 0)
# and      ( 4  3  0) 
#      A = (-2 -1  0)
#          ( 0  2,-1)
# and possibly the equality constraint  3x1 + 2x2 + x3 = 1
# or upper bound c(1.5, 1.5, 1.5).

## quadprog(C, d, A = NULL, b = NULL, Aeq = NULL, beq = NULL, lb = NULL, ub = NULL)

#install.packages("installr")
#library(installr)
#updateR()

install.packages("quadprog")
library(quadprog)
install.packages("qrmdata")
library(qrmdata)

rm(list=ls())

############# Part I. Introductions to R's quadprog library ###########################
C <- diag(1, 3)
d <- -c(0, 5, 0)
A <- matrix(c(4,3,0, -2,-1,0, 0,2,-1), 3, 3, byrow=TRUE)
b <- c(8, -2, 0)

#1
solve.QP(C, d, A, b)

#2
Aeq <- c(3, 2, 1) 
beq <- 1
solve.QP(C, d, A, b, Aeq, beq)

#3
solve.QP(C, d, A, b, lb = 0, ub = 1.5)


############################################## Part II. Example Data Preparation ###############################################
## solve.QP's arguments are:
## solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)

# Load sample data
library(quadprog)
## install.packages("qrmdata")
library(qrmdata)
require(qrmdata)
require(xts)

## The exchange rate data was obtained from OANDA (http://www.oanda.com/) on 2016-01-03
data("EUR_USD")
data("GBP_USD")

## The Brent data was obtained from Federal Reserve Economic Data (FRED) via Quandl on 2016-01-03
data("OIL_Brent")
data.1 <- na.omit(merge(EUR_USD, GBP_USD, OIL_Brent))

## Our dataset is named "R"
R <- na.omit(diff(log(data.1)) * 100)
names.R <- c("EUR.USD", "GBP.USD", "OIL.Brent")
colnames(R) <- names.R

dim(R)
head(R)

# First we compute the mean return
mean.R <- apply(R, 2, mean)
# "2" is the dimension of the dataset R, here R is a dataframe (Matrix), thus we use "2".
summary(R)
# Means are much less than medians.

# calculate the variance-covariance matrix. 
# The diagonals of this matrix are the variances, 
# so that the square root of the diagonal will yield standard deviations.
cov.R <- cov(R)

# The diagonals of this matrix are the variances, 
# so that the square root of the diagonal will yield standard deviations.
sd.R <- sqrt(diag(cov.R))
# remember these are in daily percentages

############################################## Part III. Optimization Example Starts here ######################
################ Example 1 we build the "efficient frontier" with Short selling #########################
library(quadprog)
## solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)

## set the equality constraints matrix Amat
Amat <- cbind(rep(1, 3), mean.R)
t(Amat)

## set of 300 possible target portfolio returns
mu.P <- seq(min(mean.R - 5e-04), max(mean.R + 5e-04), length = 300) 

# just to see the very first bvec. We have bvec = c(1, mu.P[i] in the following for loop.
bvec = c(1, mu.P[1])

## three other major inputs
Dmat = 2 * cov.R 
dvec = rep(0, 3)
meq=2

## the major components of the for loop
result = solve.QP(Dmat = 2 * cov.R, 
                  dvec = rep(0, 3), Amat = Amat, 
                  bvec = bvec, meq = 2)
## sigma.P[i] = sqrt(result$value)
## weights[i, ] = result$solution


## set up storage for std dev's of portfolio returns
sigma.P <- mu.P  

## storage for portfolio weights
weights <- matrix(0, nrow = 300, ncol = ncol(R))  
colnames(weights) <- names.R


## This curve (a parabola.) traces optimal combinations of risk and return. 
## For each combination there is an underlying set of weights calculated in successive optimizations, 
## one for each target u(mu)
## Next we build the "efficient frontier."
for (i in 1:length(mu.P)) {
  bvec = c(1, mu.P[i])  ## constraint vector
  result = solve.QP(Dmat = 2 * cov.R, 
                    dvec = rep(0, 3), Amat = Amat, 
                    bvec = bvec, meq = 2)
  sigma.P[i] = sqrt(result$value)
  weights[i, ] = result$solution
}

############### II. Then we plot #########################
## 1. Plot all of the portfolio combinations.
## 2. Plot the point on the graph that represents the risk-free asset.

par(mfrow = c(1, 1))
##  plot the efficient frontier (and inefficient portfolios below the
## min var portfolio)
plot(sigma.P, mu.P, type = "l", xlim = c(0, max(sd.R) * 1.1), ylim = c(0, max(mean.R) * 1.1), lty = 3, lwd = 3)

## input value of risk-free interest rate
mu.free = 1.3/253  
## show risk-free asset
points(0, mu.free, cex = 1, pch = "+")  

############### III. William Sharpe's ratio, again as before
## 1. This number is the amount of portfolio premium per unit of risk (the "price" of risk) 
## across all combinations of portfolio assets on the efficient frontier. 
## Its maximum is the best combination for the risk in terms of returns.

## 2. We figure out where (the index ind) the return to risk is along the frontier, 
## record the weights associated with this unique point in risk-return space, and

## 3. Find where (the index ind2) the minimum variance portfolio is.

## 4. Plot the "efficient frontier": the efficient frontier will extend from the minimum variance portfolio 
## (a "+" will mark the spot) up and out (in red). 
## Anything else below this line is "inefficient" in the sense you get less and less return for more and m

par(mfrow = c(1, 1))
##  plot the efficient frontier (and inefficient portfolios below the min var portfolio)
plot(sigma.P, mu.P, type = "l", xlim = c(0, 
                                         max(sd.R) * 1.1), ylim = c(min(mean.R) * 
                                                                      1.05, max(mean.R) * 1.1), lty = 3, 
     lwd = 3)  

## input value of risk-free interest rate
mu.free <- 1.3/253  

## show risk-free asset
points(0, mu.free, cex = 1.5, pch = "+")  

## input value of risk-free interest rate
mu.free <- 1.3/253

## show risk-free asset
points(0, mu.free, cex = 1.5, pch = "+")  

## compute Sharpe's ratios
sharpe = (mu.P - mu.free)/sigma.P  

## Find maximum Sharpe's ratio
ind = (sharpe == max(sharpe))  
options(digits = 3)
## show line of optimal portfolios
lines(c(0, 2), mu.free + c(0, 2) * (mu.P[ind] - 
                                      mu.free)/sigma.P[ind], lwd = 4, lty = 1, 
      col = "blue")

## show tangency portfolio
points(sigma.P[ind], mu.P[ind], cex = 4, 
       pch = "*")  

## find the minimum variance portfolio
ind2 = (sigma.P == min(sigma.P))

## show min var portfolio
points(sigma.P[ind2], mu.P[ind2], cex = 1.5, 
       pch = "+")  
ind3 = (mu.P > mu.P[ind2])

##  plot the efficient frontier
lines(sigma.P[ind3], mu.P[ind3], type = "l", 
      xlim = c(0, max(sd.R) * 1.1), ylim = c(min(mean.R) * 
                                               1.05, max(mean.R) * 1.1), lwd = 3, 
      col = "red")  
text(sd.R[1], mean.R[1], "EUR.USD", cex = 1.15)
text(sd.R[2], mean.R[2], "GBP.USD", cex = 1.15)
text(sd.R[3], mean.R[3], "OIL.Brent", 
     cex = 1.15)

## The weights for the tangency portfolio ("*")
weights[ind, ]
sum(weights[ind, ])

## Summary:
## Real world meaning: For a given notional amount in your portfolio, 
## go long (buy) 250.% of that position in euros traded against USD, 
## go short (sell) 180.7% of your aggregate position in euros traded against USD, and 
## go long 30.6% in Brent.


################ Example 2 Short selling is NOT ALLOWED #########################
## In this scenario we don't allow short positions (negative weights). 
## This means we impose the inequality constraint:
## w>=0

## Further, 
## We modify the Amat to Amat to cbind(rep(1,3),mean.R,diag(1,nrow=3)).
## We set the target return vectormu.P to seq(min(mean.R)+.0001, max(mean.R)-.0001, length=300).
## We also set the righthand-side vector bvec to c(1,mu.P[i],rep(0,3))

############### I. Here is the new setup code where we no longer allow for short positions.
library(quadprog)
## set the equality ND inequality constraints matrix
Amat <- cbind(rep(1, 3), mean.R, diag(1, nrow = 3))  
t(Amat)

bvec <- c(1, mu.P[1], rep(0, 3))

## set of 300 possible target portfolio returns
mu.P <- seq(min(mean.R) + 1e-04, max(mean.R) - 1e-04, length = 300)  

## set up storage for std dev's of portfolio returns
sigma.P <- mu.P  

## storage for portfolio weights
weights <- matrix(0, nrow = 300, ncol = 3)  

## Next we build the "efficient frontier." All of this code is as before.
for (i in 1:length(mu.P)) {
  bvec <- c(1, mu.P[i], rep(0, 3))  ## constraint vector with no short positions
  result <- solve.QP(Dmat = 2 * cov.R, 
                     dvec = rep(0, 3), Amat = Amat, 
                     bvec = bvec, meq = 2)
  sigma.P[i] <- sqrt(result$value)
  weights[i, ] <- result$solution
}

############### II. Then we plot, again the same as before
## 1. Plot all of the portfolio combinations.
## 2. Plot the point on the graph that represents the so-called risk-free (actually more like default-free) asset.

par(mfrow = c(1, 1))

## plot the efficient frontier (and inefficient portfolios
plot(sigma.P, mu.P, type = "l", xlim = c(0, 
                                         max(sd.R) * 1.1), ylim = c(min(mean.R) * 
                                                                      1.05, max(mean.R) * 1.1), lty = 3, 
     lwd = 3)  
## below the min var portfolio)

## input value of risk-free interest rate
mu.free <- 1.3/253  

## show risk-free asset
points(0, mu.free, cex = 1.5, pch = "+")  

############### III. William Sharpe's ratio, again as before
## 1. This number is the amount of portfolio premium per unit of risk (the "price" of risk) 
## across all combinations of portfolio assets on the efficient frontier. 
## Its maximum is the best combination for the risk in terms of returns.

## 2. We figure out where (the index ind) the return to risk is along the frontier, 
## record the weights associated with this unique point in risk-return space, and

## 3. Find where (the index ind2) the minimum variance portfolio is.

## 4. Plot the "efficient frontier": the efficient frontier will extend from the minimum variance portfolio 
## (a "+" will mark the spot) up and out (in red). 
## Anything else below this line is "inefficient" in the sense you get less and less return for more and m

par(mfrow = c(1, 1))
##  plot the efficient frontier (and inefficient portfolios below the min var portfolio)
plot(sigma.P, mu.P, type = "l", xlim = c(0, 
                                         max(sd.R) * 1.1), ylim = c(min(mean.R) * 
                                                                      1.05, max(mean.R) * 1.1), lty = 3, 
     lwd = 3)  

## input value of risk-free interest rate
mu.free <- 1.3/253  

## show risk-free asset
points(0, mu.free, cex = 1.5, pch = "+")  

## input value of risk-free interest rate
mu.free <- 1.3/253

## show risk-free asset
points(0, mu.free, cex = 1.5, pch = "+")  

## compute Sharpe's ratios
sharpe = (mu.P - mu.free)/sigma.P  

## Find maximum Sharpe's ratio
ind = (sharpe == max(sharpe))  
options(digits = 3)
## show line of optimal portfolios
lines(c(0, 2), mu.free + c(0, 2) * (mu.P[ind] - 
                                      mu.free)/sigma.P[ind], lwd = 4, lty = 1, 
      col = "blue")

## show tangency portfolio
points(sigma.P[ind], mu.P[ind], cex = 4, 
       pch = "*")  

## find the minimum variance portfolio
ind2 = (sigma.P == min(sigma.P))

## show min var portfolio
points(sigma.P[ind2], mu.P[ind2], cex = 1.5, 
       pch = "+")  
ind3 = (mu.P > mu.P[ind2])

##  plot the efficient frontier
lines(sigma.P[ind3], mu.P[ind3], type = "l", 
      xlim = c(0, max(sd.R) * 1.1), ylim = c(min(mean.R) * 
                                               1.05, max(mean.R) * 1.1), lwd = 3, 
      col = "red")  
text(sd.R[1], mean.R[1], "EUR.USD", cex = 1.15)
text(sd.R[2], mean.R[2], "GBP.USD", cex = 1.15)
text(sd.R[3], mean.R[3], "OIL.Brent", 
     cex = 1.15)

## Now the constraints are:
Amat
bvec

## 1. bvec changes for each of the three assets. Here we see one of them.
## 2. The short position bvec has three zeros appended to it.
## 3. The Amat constraint matrix has the identity matrix appended to it to represent:  
## wi=0
## in the formulation of the inequality constraints parsed by quantprog.
## 4. The tangency of the line from the risk-free rate to the maximum Sharpe ratio point 
## on the efficient frontier does not change.

## The weights are:
weights[ind, ]
sum(weights[ind, ])

## Summary:
## 1. Long working capital position with only a $1 million euro exposure.
## 2. No pounding sterling exposure at all.
## 3. A huge $99 million Brent exposure.