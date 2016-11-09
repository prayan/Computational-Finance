# M237G Wei Wei
###### p1
# Ch2 P38
n = 1000
a = -.7
sigma1 = 1
sigma2 = 1
rho = a/sigma1/sigma2

# covmtx = matrix( c(1, a, a, 1), ncol = 2, byrow = FALSE )

# bivariate normal distribution x
binorm = function(n, sigma1, sigma2, rho){
  z1 = rnorm(n, mean = 0, sd = 1)
  z2 = rnorm(n, mean = 0, sd = 1)
  z = matrix(c(z1, z2), nrow = 2, byrow = TRUE)
  
  L = matrix( c(sigma1, sigma2*rho, 0, sqrt(1-rho^2)*sigma2), ncol = 2, byrow = FALSE )
  
  x = matrix( NA, nrow = nrow(z), ncol = ncol(z))
  for( i in 1:ncol(x) ){
    x[,i] = L %*% z[,i]
  }
  x
}

x = binorm(n, sigma1, sigma2, rho)

# verify
mu_x1 = mean(x[1,])
mu_x2 = mean(x[2,])
# sample covariance
numerator = sum( (x[1,]-mu_x1) * (x[2,]- mu_x2) , na.rm = FALSE) / (n-1)
# sqrt(sample variance)
denominator = sqrt(var(x[1,]) * var(x[2,]))
rho_sim = numerator / denominator


###### p2
xy = binorm(n, 1, 1, 0.6)
E = (x[1,]^3 + sin(x[2,]) + x[1,]^2 * x[2,])/n
E = max (0, E)

###### p3
Z = rnorm(n, 0, 1)
w5 = sqrt(5) * Z
Ea1 = mean(w5^2 + sin(w5))
var(w5^2 + sin(w5))

expectation = function(t, Z){
  W_t = sqrt(t) * Z
  E = mean(exp(t/2) * cos(W_t))
  return (E)
}
Ea2 = expectation(t = 0.5,Z)
Ea3 = expectation(t = 3.2,Z)
Ea4 = expectation(t = 6.5,Z)

#### 3.(b) Variance Reduction
z1 = rnorm(n, mean = 0, sd = 1)
z2 = -z1
w5b = sqrt(5) * z1    # W(5)+
w5b_ = sqrt(5) * z2   # W(5)_
Eb1 = (mean(w5b^2 + sin(w5b)) + mean(w5b_^2 + sin(w5b_)) )/2
Eb2 = (expectation(t = 0.5, z1) + expectation(t = 0.5, z2))/2
Eb3 = (expectation(t = 3.2, z1) + expectation(t = 3.2, z2))/2
Eb4 = (expectation(t = 6.5, z1) + expectation(t = 6.5, z2))/2

###### p4
r = .04
sigma = .2
s0 = 88
maturity = 5
K = 100
n = 1000

# European call price single simulation
eurocall = function(s0, r, sigma, maturity, K, W) {
  price = max( s0 * exp( (r-0.5*sigma^2)*maturity + sigma*W ) - K, 0 )
  return (price)
}

# European call Black-Scholes price
myBSPrice <- function( s0, K, r, sigma, maturity ) {
  d1 = (log(s0/K)+(r+sigma^2/2)*maturity)/(sigma*sqrt(maturity))
  d2 = (log(s0/K)+(r-sigma^2/2)*maturity)/(sigma*sqrt(maturity))
  
  # pnorm is CDF
  call_price <- s0*pnorm(d1)-exp(-r*maturity)*K*pnorm(d2)
  call_price
}

Z = rnorm(n,mean = 0, sd = 1)
# Brownian Motion
W = sqrt(maturity) * Z

# call price by Monte Carlo simulation
payoff = array(NaN, dim = n)
for( i in 1:length(W) ){
  payoff[i] = eurocall(s0, r, sigma, maturity, K, W[i])
}
Ca1 = exp(-r*maturity) * mean(payoff)

# call price by BS formula
Cb1 = myBSPrice(s0,K,r,sigma,maturity)

# call price from Monte Carlo sim with Variance Reduction
payoff_pos = array(0, dim = n)
payoff_neg = array(0, dim = n)
z1 = rnorm(n,mean = 0, sd = 1)
z2 = -z1
for( i in 1:length(z1)){
  payoff_pos[i] = eurocall(s0, r, sigma, maturity, K, sqrt(maturity) * z1[i])
}
for( i in 1:length(z2)){
  payoff_neg[i] = eurocall(s0, r, sigma, maturity, K, sqrt(maturity) * z2[i])
}
exp(-r*maturity) * mean(payoff_pos)
exp(-r*maturity) * mean(payoff_neg)
Cb2 = exp(-r*maturity) * (mean(payoff_pos) + mean(payoff_neg))/2

###### p5
# Geometric Brownian Motion
myGBM = function(r, sigma, s0, t){
  W = sqrt(t) * rnorm(1000, 0, 1)
  st = s0 * mean(exp( (r-0.5*sigma^2)*t + sigma*W ))
  st
}
r = 0.04
sigma = 0.18
s0 = 88
ESn = array(0, dim = 10)
for( i in 1: length(ESn)){
  ESn[i] = myGBM(r, sigma, s0, t = i)
}
plot(ESn, type = "l", col = "red")

# 5(b)
t = (1:1000)/100
ESn_b = matrix(NA, nrow = 1000, ncol = 6)
for( i in 1:6 ){
  for( j in 1:1000){
    ESn_b[j,i] = myGBM(r, sigma, s0, t[j])
  }
}

for( i in 1:6 ){
  par( new = TRUE )
  plot(x = t, y = ESn_b[,i], type = "l")
}

# 5(d)
sigma2 = 0.35
ESn_5d1 = array(0, dim = 10)
for( i in 1: length(ESn_5d1)){
  ESn_5d1[i] = myGBM(r, sigma, s0, t = i)
}
plot(ESn, type = "l", col = "red")
par( new = TRUE)
plot(ESn_5d1, type = "l", col = "green", ylab = "")

ESn_5d2 = matrix(NA, nrow = 1000, ncol = 6)
for( i in 1:6 ){
  for( j in 1:1000){
    ESn_5d2[j,i] = myGBM(r, sigma, s0, t[j])
  }
}

for( i in 1:6 ){
  par( new = TRUE )
  plot(x = t, y = ESn_5d2[,i], type = "l", ylab = "")
}

###### p6
# (a) Euler's discretization
x = (1:10000)/10000
Ia = sum(4 * sqrt(1-x^2) * (0.0001))
# (b) Monte Carlo
x_6b = runif(10000, 0, 1)     # x~U[0,1]
Ib = 4 * mean(sqrt(1-x_6b^2))
var(sqrt(1-x_6b^2))

# (c) importance sampling
# t(x)
myT = function(x){
  a = 0.74
  if(x>=0 && x<=1){
    t = (1-a*x^2) / (1 - a/3)
  } else{
    t = 0
  }
  return (t)
}
# define g(x) max(t(x)) = t(0) = 1.327434
myT(0)

U = runif(1000, 0, 1)
# if t(x) / g(x) >= U, accept t(x)
Y = (((1-a*U^2)/(1-a/3)/1.327434) >= U) * U

Values = array(0, dim = length(Y))
for( i in length(Y)){
  Values[i] = 4 * sqrt(1 - Y[i]^2) / (1-a * Y[i]^2) * (1-a/3)
}
mean(Values)
Ic = mean(4 * sqrt(1 - Y^2) / (1-a * Y^2) * (1-a/3))
var(sqrt(1 - Y^2) / (1-a * Y^2) * (1-a/3))
