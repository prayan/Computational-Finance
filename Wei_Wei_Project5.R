# Wei Wei MGMT237G Project 5
# problem 1 ####################################
sigma = 0.2
rf = 0.06
paths = 100000 # in computation, default number of paths is 60 just to check results.
X = 40
mT = c(0.5,1,2)
s0 = c(36,40,44)

### common function definition ###

# Returning the first >0 value in each row of a matrix
firstValueRow = function(x) {
  cumSumMat = matrix(NA, nrow=dim(x)[1], ncol=dim(x)[2])
  for(i in 1:(dim(x)[1])) {
    cumSumMat[i,] = cumsum(x[i,])
  }
  cumSumMat2 = cbind(matrix(0, nrow=dim(x)[1], ncol=1), cumSumMat[,-(dim(x)[2])])
  ResultMat = matrix(NA, nrow=dim(x)[1], ncol=dim(x)[2])
  for(i in 1:dim(x)[2]) {
    ResultMat[,i] = ifelse(cumSumMat2[,i]>0, 0, x[,i])
  }
  return(ResultMat)
}

summary.AmerPut = function(object){
  cat("\nAmerican Put Option\nMethod: Simple Least Squares Monte Carlo\n")
  mat1 = matrix(object)
  rownames(mat1) = c("Option price", "Spot price", "Strike", "Volatility", "Number of paths", "Number of time-steps", "Interest rate", "Dividend rate", "Maturity time")
  colnames(mat1) = c(" ")
  print(mat1)
}

### 1(a) Laguerre polynomials for k = 2,3,4
### function definition ###
AmerPutLSM_Lag = function(Spot=1, sigma=0.2, n=1000, m=365, Strike=1.1, r=0.06, dr=0.0, mT=1, k=4) { # n paths, m time-steps, k terms
  GBM = matrix(NA, nrow=n, ncol=m)
  for(i in 1:n) {
    GBM[i,] = Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X = ifelse(GBM<Strike,GBM,NA)
  CFL = matrix(pmax(0,Strike-GBM), nrow=n, ncol=m)
  # Laguerre polynomials
  Xsh = exp(-X[,-m]/2)
  if( k == 3 ){
    X2sh = Xsh*(1-X[,-m])
  }
  if( k == 4 ){
    X2sh = Xsh*(1-X[,-m])
    X3sh = Xsh * (1-2*X[,-m]+0.5*X[,-m]^2)
  }
  Y1 = CFL*exp(-1*r*(mT/m))
  Y2 = cbind((matrix(NA, nrow=n, ncol=m-1)), Y1[,m])
  CV = matrix(NA, nrow=n, ncol=m-1)
  try(for(i in (m-1):1) {
    if( k == 2 ){
      reg1 = lm(Y2[,i+1]~Xsh[,i])
      CV[,i] = (matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])
    }
    else if( k == 3 ){
      reg1 = lm(Y2[,i+1]~Xsh[,i]+X2sh[,i])
      CV[,i] = (matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])
    }
    else if( k == 4 ){
      reg1 = lm(Y2[,i+1]~Xsh[,i]+X2sh[,i]+X3sh[,i])
      CV[,i] = (matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])+((matrix(reg1$coefficients)[4,1])*X3sh[,i])
    }
    CV[,i] = (ifelse(is.na(CV[,i]),0,CV[,i]))
    Y2[,i] = ifelse(CFL[,i]>CV[,i], Y1[,i], Y2[,i+1]*exp(-1*r*(mT/m)))
  }
  , silent = TRUE)
  CV = ifelse(is.na(CV),0,CV)
  CVp = cbind(CV, (matrix(0, nrow=n, ncol=1)))
  POF = ifelse(CVp>CFL,0,CFL)
  FPOF = firstValueRow(POF)
  dFPOF = matrix(NA, nrow=n, ncol=m)
  for(i in 1:m) {
    dFPOF[,i] = FPOF[,i]*exp(-1*mT/m*r*i)
  }
  PRICE = mean(rowSums(dFPOF))
  res =  list(price=(PRICE), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res) = "AmerPut"
  return(res)
}
### computation ###
# k=2
result_ak2 = matrix(NA, nrow = length(s0), ncol = length(mT))
for( i in 1:length(s0) ){
  for( j in 1:length(mT) ){
    amerput = AmerPutLSM_Lag(Spot = s0[i], sigma, n = 60, Strike = X, r = rf, mT = mT[i], k=2 )
    summary.AmerPut(amerput)
    result_ak2[i,j] = amerput[[1]]
  }
}
# k = 3
result_ak3 = matrix(NA, nrow = length(s0), ncol = length(mT))
for( i in 1:length(s0) ){
  for( j in 1:length(mT) ){
    amerput = AmerPutLSM_Lag(Spot = s0[i], sigma, n = 60, Strike = X, r = rf, mT = mT[i], k=3 )
    summary.AmerPut(amerput)
    result_ak3[i,j] = amerput[[1]]
  }
}
# k = 4
result_ak4 = matrix(NA, nrow = length(s0), ncol = length(mT))
for( i in 1:length(s0) ){
  for( j in 1:length(mT) ){
    amerput = AmerPutLSM(Spot = s0[i], sigma, n = 60, Strike = X, r = rf, mT = mT[i], k=4 )
    summary.AmerPut(amerput)
    result_ak4[i,j] = amerput[[1]]
  }
}

########################################
# 1(b) Hermite polynomials for k = 2,3,4
### function definition ###
AmerPutLSM = function(Spot=1, sigma=0.2, n=1000, m=365, Strike=1.1, r=0.06, dr=0.0, mT=1, k=4) { # n paths, m time-steps, k terms
  GBM = matrix(NA, nrow=n, ncol=m)
  for(i in 1:n) {
    GBM[i,] = Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X = ifelse(GBM<Strike,GBM,NA)
  CFL = matrix(pmax(0,Strike-GBM), nrow=n, ncol=m)
  # Hermite polynomials
  Xsh = 2*X[,-m]
  if( k == 3 ){
    X2sh = Xsh*Xsh - 2
  }
  if( k == 4 ){
    X2sh = Xsh*Xsh - 2
    X3sh = Xsh*Xsh*Xsh - 12 * X[,-m]
  }
  Y1 = CFL*exp(-1*r*(mT/m))
  Y2 = cbind((matrix(NA, nrow=n, ncol=m-1)), Y1[,m])
  CV = matrix(NA, nrow=n, ncol=m-1)
  try(for(i in (m-1):1) {
    if( k == 2 ){
      reg1 = lm(Y2[,i+1]~Xsh[,i])
      CV[,i] = (matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])
    }
    else if( k == 3 ){
      reg1 = lm(Y2[,i+1]~Xsh[,i]+X2sh[,i])
      CV[,i] = (matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])
    }
    else if( k == 4 ){
      reg1 = lm(Y2[,i+1]~Xsh[,i]+X2sh[,i]+X3sh[,i])
      CV[,i] = (matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])+((matrix(reg1$coefficients)[4,1])*X3sh[,i])
    }
    CV[,i] = (ifelse(is.na(CV[,i]),0,CV[,i]))
    Y2[,i] = ifelse(CFL[,i]>CV[,i], Y1[,i], Y2[,i+1]*exp(-1*r*(mT/m)))
  }
  , silent = TRUE)
  CV = ifelse(is.na(CV),0,CV)
  CVp = cbind(CV, (matrix(0, nrow=n, ncol=1)))
  POF = ifelse(CVp>CFL,0,CFL)
  FPOF = firstValueRow(POF)
  dFPOF = matrix(NA, nrow=n, ncol=m)
  for(i in 1:m) {
    dFPOF[,i] = FPOF[,i]*exp(-1*mT/m*r*i)
  }
  PRICE = mean(rowSums(dFPOF))
  res =  list(price=(PRICE), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res) = "AmerPut"
  return(res)
}

### computation ###
# k=2
result_bk2 = matrix(NA, nrow = length(s0), ncol = length(mT))
for( i in 1:length(s0) ){
  for( j in 1:length(mT) ){
    amerput = AmerPutLSM_Lag(Spot = s0[i], sigma, n = 60, Strike = X, r = rf, mT = mT[i], k=2 )
    summary.AmerPut(amerput)
    result_bk2[i,j] = amerput[[1]]
  }
}
# k = 3
result_bk3 = matrix(NA, nrow = length(s0), ncol = length(mT))
for( i in 1:length(s0) ){
  for( j in 1:length(mT) ){
    amerput = AmerPutLSM_Lag(Spot = s0[i], sigma, n = 60, Strike = X, r = rf, mT = mT[i], k=3 )
    summary.AmerPut(amerput)
    result_bk3[i,j] = amerput[[1]]
  }
}
# k = 4
result_bk4 = matrix(NA, nrow = length(s0), ncol = length(mT))
for( i in 1:length(s0) ){
  for( j in 1:length(mT) ){
    amerput = AmerPutLSM(Spot = s0[i], sigma, n = 60, Strike = X, r = rf, mT = mT[i], k=4 )
    summary.AmerPut(amerput)
    result_bk4[i,j] = amerput[[1]]
  }
}

#####################################
# 1(c) simple monomials for k = 2,3,4
### function definition ###
AmerPutLSM = function(Spot=1, sigma=0.2, n=1000, m=365, Strike=1.1, r=0.06, dr=0.0, mT=1, k=4) { # n paths, m time-steps, k terms
  GBM = matrix(NA, nrow=n, ncol=m)
  for(i in 1:n) {
    GBM[i,] = Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  X = ifelse(GBM<Strike,GBM,NA)
  CFL = matrix(pmax(0,Strike-GBM), nrow=n, ncol=m)
  # simple monomials
  Xsh = X[,-m]
  if( k == 3 ){
    X2sh = Xsh*Xsh
  }
  if( k == 4 ){
    X2sh = Xsh*Xsh
    X3sh = Xsh*Xsh*Xsh
  }
  Y1 = CFL*exp(-1*r*(mT/m))
  Y2 = cbind((matrix(NA, nrow=n, ncol=m-1)), Y1[,m])
  CV = matrix(NA, nrow=n, ncol=m-1)
  try(for(i in (m-1):1) {
    if( k == 2 ){
      reg1 = lm(Y2[,i+1]~Xsh[,i])
      CV[,i] = (matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])
    }
    else if( k == 3 ){
      reg1 = lm(Y2[,i+1]~Xsh[,i]+X2sh[,i])
      CV[,i] = (matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])
    }
    else if( k == 4 ){
      reg1 = lm(Y2[,i+1]~Xsh[,i]+X2sh[,i]+X3sh[,i])
      CV[,i] = (matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])+((matrix(reg1$coefficients)[4,1])*X3sh[,i])
    }
    CV[,i] = (ifelse(is.na(CV[,i]),0,CV[,i]))
    Y2[,i] = ifelse(CFL[,i]>CV[,i], Y1[,i], Y2[,i+1]*exp(-1*r*(mT/m)))
  }
  , silent = TRUE)
  CV = ifelse(is.na(CV),0,CV)
  CVp = cbind(CV, (matrix(0, nrow=n, ncol=1)))
  POF = ifelse(CVp>CFL,0,CFL)
  FPOF = firstValueRow(POF)
  dFPOF = matrix(NA, nrow=n, ncol=m)
  for(i in 1:m) {
    dFPOF[,i] = FPOF[,i]*exp(-1*mT/m*r*i)
  }
  PRICE = mean(rowSums(dFPOF))
  res =  list(price=(PRICE), Spot, Strike, sigma, n, m, r, dr, mT)
  class(res) = "AmerPut"
  return(res)
}

### computation ###
# k=2
result3 = matrix(NA, nrow = length(s0), ncol = length(mT))
for( i in 1:length(s0) ){
  for( j in 1:length(mT) ){
    amerput = AmerPutLSM(Spot = s0[i], sigma, n = 60, Strike = X, r = rf, mT = mT[i], k=2 )
    summary.AmerPut(amerput)
    result3[i,j] = amerput[[1]]
  }
}
# k = 3
result3 = matrix(NA, nrow = length(s0), ncol = length(mT))
for( i in 1:length(s0) ){
  for( j in 1:length(mT) ){
    amerput = AmerPutLSM(Spot = s0[i], sigma, n = 60, Strike = X, r = rf, mT = mT[i], k=3 )
    summary.AmerPut(amerput)
    result3[i,j] = amerput[[1]]
  }
}
# k = 4
result3 = matrix(NA, nrow = length(s0), ncol = length(mT))
for( i in 1:length(s0) ){
  for( j in 1:length(mT) ){
    amerput = AmerPutLSM(Spot = s0[i], sigma, n = 60, Strike = X, r = rf, mT = mT[i], k=4 )
    summary.AmerPut(amerput)
    result3[i,j] = amerput[[1]]
  }
}

# problem 2 ####################################
# 2(a) forward start European put option pricing
s0 = 65
sigma = 0.2
rf = 0.06
t = 0.2
mT = 1 # maturity
n = 1000 # number of stock price paths
m = 60 # number of time steps from 0 to t

### function definition ###

# brownian motion path
fastGBM  = function(Spot=1, sigma=0.2, n=1000, m=365, r=0.06, dr=0.0, mT=1) {
  GBM = matrix(NA, nrow=n, ncol=m)
  for(i in 1:n) {
    GBM[i,] = Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
  }
  return(GBM)
}

# European put black-scholes price
BSPut_euro = function(s0, K, r, sigma, mT){
  d1 = (log(s0/K) + (r+0.5*sigma^2)*mT) / (sigma*sqrt(mT))
  d2 = d1 - sigma*sqrt(mT)
  K*exp(-r*mT)*pnorm(-d2) - s0*pnorm(-d1)
}

### computation ###

# simulate 1000 paths of stock price Brownian Motion from time 0 to t, number of stops is 60
st = fastGBM(Spot = s0, sigma, n, m, r = rf, dr=0, mT = t)
# at strike determination time t, compute standard BS put price for each stock price path
put_euro_t = BSPut_euro(st[,m], st[,m], rf, sigma, mT-t)
# take expectation of put price at time t and discount to time 0
fsput_euro = exp(-rf*t) * mean(put_euro_t)
fsput_euro

################################################
# 2(b) forward-start American put option pricing

### function definition ###

# American put
americanPut = function(s0,u,d,K,steps,r,maturity){
  step = maturity/steps
  R = exp(r*step)
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(s0), ncol = n+1) # option price, each row is an observation
  for( row in 1:length(s0) ){
    for( i in 0:n ){
      P[row,i+1] = max(0, K - s0[row] * u^(n-i) * d^i)
    }
    
    for( j in n:1 ){
      for( k in 1:j ){
        P[row,k] = max( (p*P[row,k] + (1-p)*P[row,k+1])/R, K - s0 * u^(j-k+1) * d^(k-1) )
      }
    }
  }
  return ( P[,1] )
}

### computation ###

# use the stock price paths simulated in 2(a): st
# at strike determination time, compute american put price
put_amer_t = americanPut(st[m], u=1.5, d = 0.5, st[m], steps = 10000, r = rf, maturity = mT-t)
# take expectation of put price at time t and discount to time 0
fsput_amer = exp(-rf*t) * mean(put_amer_t)
fsput_amer
