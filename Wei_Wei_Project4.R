# 237G Project 4 Wei Wei

#### problem 1
r = .05
sigma = .24
s0 = 32
X = 30
maturity = 0.5
n = c(10,20,40,80,100,200,500)

binomialcall = function (u, d, p, s0, X, n, r, step) {
  R = exp(r * step)
  # c[j,i] represents the call value at time j if the stock price makes i downward movements out of a total of j movements.
  C = matrix(data = NA, nrow = n+1, ncol = n+1)
  for( i in 0:n ){
    C[n+1,i+1] = max(0, s0 * u^(n-i) * d^(i) - X)
  }
  for( j in n:1 ){
    for ( k in 1:j ){
      C[j,k] = (p * C[j+1,k] + (1-p) * C[j+1,k+1]) / R
    }
  }
  return( C[1,1] )
}

# (a)
step = maturity / n
c = 0.5 * ( exp(-r*step) + exp((r+sigma^2)*step))
d = c - sqrt(c^2-1)
u = 1/d
p = (exp(r*step) - d) / (u-d)
price_1a = array(data = NA, dim = length(n))
for (i in 1:length(n) ){
  price_1a[i] = binomialcall(u[i], d[i], p[i], s0, X, n[i], r, step[i])
}
# (b)
u = exp(r*step) * (1 + sqrt(exp(sigma^2*step)-1))
d = exp(r*step) * (1 - sqrt(exp(sigma^2*step)-1))
p = 1/2
price_1b = array(data = NA, dim = length(n))
for (i in 1:length(n) ){
  price_1b[i] = binomialcall(u[i], d[i], p, s0, X, n[i], r, step[i])
}
# (c)
u = exp((r-0.5*sigma^2)*step + sigma*sqrt(step))
d = exp((r-0.5*sigma^2)*step - sigma*sqrt(step))
p = 1/2
price_1c = array(data = NA, dim = length(n))
for (i in 1:length(n) ){
  price_1c[i] = binomialcall(u[i], d[i], p, s0, X, n[i], r, step[i])
}
# (d)
u = exp(sigma * sqrt(step))
d = exp( - sigma * sqrt(step))
p = 0.5 + 0.5*( (r-0.5*sigma^2) * sqrt(step) / sigma )
price_1d = array(data = NA, dim = length(n))
for (i in 1:length(n) ){
  price_1d[i] = binomialcall(u[i], d[i], p[i], s0, X, n[i], r, step[i])
}
par(mfrow=c(2,2))
plot(x = n, y = price_1a, xlab = "n", col = "red")
plot(x = n, y = price_1b, xlab = "n", col = "red")
plot(x = n, y = price_1c, xlab = "n", col = "red")
plot(x = n, y = price_1d, xlab = "n", col = "red")
par(mfrow=c(1,1))

#### problem 2
r = 0.02
s0 = 696.95 # GOOG stock price on April 27 2016
sigma = .24
K = floor(1.1*s0) # strike price 766
maturity = 260/360 # assume time to expiration date Jan 12 2017 is 260 days
# use the binomial call method in problem 1(c)
n = 20
step = maturity/n # divide the time interval into 100 periods.
u = exp((r-0.5*sigma^2)*step + sigma*sqrt(step))
d = exp((r-0.5*sigma^2)*step - sigma*sqrt(step))
p = 1/2
price_p2a = binomialcall(u,d,p, s0, K, n, r, step) # price = 34.12, based on sigma = .24
price_p2a
# (b)
# yahoo finance GOOG call(K = 765) price is 30.00
sigma = .2207
u = exp((r-0.5*sigma^2)*step + sigma*sqrt(step))
d = exp((r-0.5*sigma^2)*step - sigma*sqrt(step))
price_p2b = binomialcall(u,d,p, s0, K, n, r, step)
price_p2b

#### problem 3
K = 50
r = 0.03
sigma = 0.2
maturity = 0.3846 # 20 weeks
mu = 0.25
# (i)
DeltaCall = function(s0,u,d,K,n,r,maturity){
  step = maturity/n
  R = exp(r*step)
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(s0), ncol = n+1) # call option price
  for( row in 1:length(s0) ){
  for( i in 0:n ){
    P[row,i+1] = max(0, s0[row] * u^(n-i) * d^i - K )
  }
  for( j in (n:2) ){
    for( k in 1:j ){
      P[row,k] = ( p*P[row,k] + (1-p) * P[row,k+1] )/R
    }
  }
  }
  return ( (P[,1] - P[,2]) / (s0*u - s0*d) )
}

s0 = seq(20,80,2)
deltas = DeltaCall(s0, mu, 1/mu, K, 100, r, maturity) # 100 periods
plot(deltas, type = "l",main = "Delta(call)", ylab = "Delta of call option", xlab = "s0 from $20 to $80")

# (ii)
s0 = 49
maturity = seq(0,0.3846, 0.01)
step = 0.01
n = 100
DeltaCall_T =  function(s0,u,d,K,n,r,maturity){
  step = maturity / n
  R = exp(r*step)
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(maturity), ncol = n+1) # call option price, each column is a period
  for( i in 0:n ){
    P[,i+1] = max(0, s0 * u^(n-i) * d^i - K )
  }
  for( j in (n:2) ){
    for( k in 1:j ){
      P[,k] = ( p*P[,k] + (1-p) * P[,k+1] )/R
    }
  }
  return ( (P[,1] - P[,2]) / (s0*u - s0*d) )
}
deltas_T = DeltaCall_T(s0, mu, 1/mu, K, n, r, maturity)
plot(deltas_T, type = "l", xlab = "Time to expiration", ylab = "Delta")

# (iii)
Theta  = function(s0,u,d,K,n,r,maturity){
  step = maturity/n
  R = exp(r*step)
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(s0), ncol = n+1) # call option price
  for( row in 1:length(s0) ){
  for( i in 0:n ){
    P[row,i+1] = max(0, s0[row] * u^(n-i) * d^i - K )
  }
  }
  for( j in (n:2) ){
    for( k in 1:j ){
      P[,k] = ( p*P[,k] + (1-p) * P[,k+1] )/R
    }
  }
  return ( -(P[,1] - P[,2])/(step) )
}
s0 = seq(20,80,2)
maturity = 0.3846 # 20 weeks
n = 100
thetas = Theta(s0,u,d,K,n,r,maturity)
plot(thetas, main = "Theta(call)", xlab = "stock price")

# (iv)
Gamma = function(s0,u,d,K,n,r,maturity){
  step = maturity/n
  R = exp(r*step)
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(s0), ncol = n+1) # call option price
  for( row in 1:length(s0) ){
    for( i in 0:n ){
      P[row,i+1] = max(0, s0[row] * u^(n-i) * d^i - K )
    }
    for( j in (n:2) ){
      for( k in 1:j ){
        P[row,k] = ( p*P[row,k] + (1-p) * P[row,k+1] )/R
      }
    }
  }
  return ( (P[,3]+P[,1] - 2*P[,2])/(step*step) )
}
gammas = Gamma(s0,u,d,K,n,r,maturity)
plot(gammas, main = "Gamma(call)", xlab = "stock price", xlim = range(0:120))

# (v) Vega must be computed from running the binomial method twice
Vega = function(s0, sigma ,K,n,r,maturity){
  step = maturity/n
  R = exp(r*step)
  
  u = exp(sigma*sqrt(step))
  d = exp(- sigma*sqrt(step))
  
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(s0), ncol = n+1) # call option price
  for( row in 1:length(s0) ){
  for( i in 0:n ){
    P[row,i+1] = max(0, s0[row] * u^(n-i) * d^i - K )
  }
  for( j in (n:2) ){
    for( k in 1:j ){
      P[row,k] = ( p*P[row,k] + (1-p) * P[row,k+1] )/R
    }
  }
  }
  c1 = P[,1]
  
  dSigma = sigma/n
  sigma = sigma * (1+1/n)
  u = exp(sigma*sqrt(step))
  d = exp(- sigma*sqrt(step))
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(s0), ncol = n+1) # call option price
  for( row in 1:length(s0) ){
    for( i in 0:n ){
      P[row,i+1] = max(0, s0[row] * u^(n-i) * d^i - K )
    }
    for( j in (n:2) ){
      for( k in 1:j ){
        P[row,k] = ( p*P[row,k] + (1-p) * P[row,k+1] )/R
      }
    }
  }
  c2 = P[,1]
  
  Vega = (c2 - c1)/dSigma
  return(Vega)
}
vegas = Vega(s0, sigma,K,n,r,maturity)
plot(vegas, main = "Vega(call)", xlab = "stock price")

# (vi) Rho must be computed from running the binomial method twice
Rho  = function(s0,u,d,K,n,r,maturity){
  step = maturity/n
  R = exp(r*step)
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(s0), ncol = n+1) # call option price
  for( row in 1:length(s0) ){
    for( i in 0:n ){
      P[row,i+1] = max(0, s0[row] * u^(n-i) * d^i - K )
    }
    for( j in (n:2) ){
      for( k in 1:j ){
        P[row,k] = ( p*P[row,k] + (1-p) * P[row,k+1] )/R
      }
    }
  }
  p1 = P[,1]
  
  dr = r/n
  r = r*(1+1/n)
  R = exp(r*step)
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(s0), ncol = n+1) # call option price
  for( row in 1:length(s0) ){
    for( i in 0:n ){
      P[row,i+1] = max(0, s0[row] * u^(n-i) * d^i - K )
    }
    for( j in (n:2) ){
      for( k in 1:j ){
        P[row,k] = ( p*P[row,k] + (1-p) * P[row,k+1] )/R
      }
    }
  }
  p2 = P[,1]
  
  return ( (p2 - p1)/dr )
}
rhos = Rho(s0,u,d,K,n,r,maturity)
plot(rhos, main = "Rho(call)", xlab = "stock price")

#### problem 4
r = 0.05
sigma = 0.3
X = 100
maturity = 1
n = 3 # assume a 3-period binomial tree
s0 = seq(80,120,4)
europeanPut = function(s0,u,d,K,n,r,maturity){
  step = maturity/n
  R = exp(r*step)
  p = (R - d)/(u-d) # risk-neutral probability
  P = matrix(data = NA, nrow = length(s0), ncol = n+1) # option price, each row is an observation
  for( row in 1:length(s0) ){
    for( i in 0:n ){
      P[row,i+1] = max(0, K - s0[row] * u^(n-i) * d^i)
    }
  }
  for( j in n:2 ){
    for( k in 1:j ){
      P[,k] = (p*P[,k] + (1-p)*P[,k+1])/R
    }
  }
  return ( P[,1] )
}

americanPut = function(s0,u,d,K,n,r,maturity){
  step = maturity/n
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

price_p4_euro = europeanPut( s0,u=1.13,d=0.8,X,n,r,maturity)
price_p4_amer = americanPut( s0,u=1.13,d=0.8,X,n,r,maturity)
plot(price_p4_euro, col = "red")
par(new = T)
plot(price_p4_amer, col = "blue", axis = FALSE, ylab = "", xlab = "")
legend("topright", inset=.05,
       c("European","American"), fill=c("red", "blue"), horiz=TRUE)

#### problem 5
r = 0.05
sigma = 0.24
s0 = 32
X = 30
maturity = 0.5 # 6 months
n = c(10, 15, 20, 40, 70, 80, 100, 200, 500)
TrinomialEuroCall = function(s0, Pu, Pd, X, n, maturity, r){
  Pm = 1 - Pu - Pd
  R = exp(r*maturity)
  C = matrix(data = NA, nrow = length(n), ncol = 2*n[length(n)]+1)
  for( row in 1:length(n) ){
    ncol = 2*n[row]+1
    for( i in 1:ncol ){
      C[row,i] = max(0, s0*u^(n[row]-i+1)*d^(i-1) - X)
    }
    for( j in (n[row]-1):0 ){
      for( k in 1:(2*j)){
        C[row,k+1] = Pu[row]*C[row,k+1] + Pm[row]*C[row,k+2] + Pd[row]*C[row,k+3]
      }
    }
  }
  C[,1]=C[,1]/R
  return(C[,1])
}
# (a)
dt = maturity/n
d = exp(-sigma*sqrt(3*dt))
u = 1/d
Pd = (r*dt*(1-u) + (r*dt)^2 + sigma^2*dt) / (u-d) / (1-d)
Pu = (r*dt*(1-d) + (r*dt)^2 + sigma^2*dt)/(u-d)/(u-1)
price_p5a = TrinomialEuroCall(s0, Pu, Pd, X, n, maturity, r)

# (b)
deltaXu = sigma*sqrt(3*dt)
deltaXd = -sigma*sqrt(3*dt)
Pd = 0.5*( (sigma^2*dt + (r-0.5*sigma^2)^2*dt^2) / (dt*deltaXu^2) - (r - sigma^2/2)*dt / dt / deltaXu)
Pu = 0.5*( (sigma^2*dt + (r-0.5*sigma^2)^2*dt^2) / (dt*deltaXu^2) + (r - sigma^2/2)*dt / dt / deltaXu)
price_p5b = TrinomialEuroCall(s0, Pu, Pd, X, n, maturity, r)

plot(price_p5a, ylab="Call price", col = "red")
legend("topright", inset=.05,
       c("a","b"), fill=c("red", "blue"), horiz=TRUE)
par(new = TRUE)
plot(price_p5b, ylab = "", axes = FALSE, col = "blue")

#### problem 6
GetHalton = function(HowMany, Base){
  seq = array(data = NA, dim = HowMany )
  NumBits = 1 + ceiling( log(HowMany)/log(Base) ) # number of points
  VecBase = Base^(-(1:NumBits)) # base^(-i) vector
  WorkVec = array(0, dim = NumBits)
  for( i in 1:HowMany ){
    j = 1
    ok = 0
    while( ok == 0 ){
      WorkVec[j] = WorkVec[j]+1
      if(WorkVec[j]<Base){
        ok = 1
      }
      else{
        WorkVec[j] = 0
        j = j+1
      }
    }
    seq[i] = sum(WorkVec * VecBase)
  }
  return(seq)
}

haltonEuropeanCall = function(s0, K, maturity, r, sigma, N, b1, b2){
  h1 = GetHalton(N, b1)
  h2 = GetHalton(N, b2)
  z1 = sqrt(-2*log(h1)) * cos(2*3.14*h2)
  z2 = sqrt(-2*log(h1)) * sin(2*3.14*h2)
  WT1 = sqrt(maturity) * z1
  WT2 = sqrt(maturity) * z2
  c1 = mean(max(0, s0 * exp((r-0.5*sigma^2)*maturity + sigma*WT1) - K )) * exp(-r*maturity)
  c2 = mean(max(0, s0 * exp((r-0.5*sigma^2)*maturity + sigma*WT2) - K )) * exp(-r*maturity)
  c = (c1+c2)/2
}

s0 = 80
K = 100
r = 0.05
maturity = 1
sigma = 0.24
N = 1000
b1 = 2
b2 = 3
C = haltonEuropeanCall(s0, K, maturity, r, sigma, N, b1, b2)