# M237G Project 8 Fixed Income pricing
# Wei Wei

nSims = 1000 # 100 Monte Carlo simulations
#############################
# Problem 1 Vasicek model, Monte Carlo, bond pricing
r0 = 0.05
sigma = 0.1
kappa = 0.82
rhat = 0.05

#### FUNCTION DEFINITION ####
Vasicek = function(t, mT, r0, k, rhat, sigma){
  # each time step is a day, assume 360 days a year
  dt = 1/360
  nstep = (mT-t) / dt
  
  r = matrix(NA, nrow = nSims, ncol = nstep) # 100 Monte Carlo simulations
  r[,1] = r0
  for( i in 2:nstep ){
    dw = rnorm(nSims)
    r[,i] = r[,i-1] + k * (rhat - r[,i-1]) * dt + sigma * sqrt(dt) * dw
  }
  integration = apply(r, 1, sum) * dt
  return (integration)
}

Vasicek.rt = function(t, mT, r0, k, rhat, sigma){
  # each time step is a day, assume 360 days a year
  dt = 1/360
  nstep = (mT-t) / dt
  
  r = matrix(NA, nrow = nSims, ncol = nstep) # nSims of Monte Carlo simulations
  r[,1] = r0
  for( i in 2:nstep ){
    dw = rnorm(nSims)
    r[,i] = r[,i-1] + k * (rhat - r[,i-1]) * dt + sigma * sqrt(dt) * dw
  }
  
  return (r[,nstep])
}

bondprice = function(t, mT, r0, k, rhat, sigma, coupon, style = c("zero", "coupon")){
  if(style=="zero"){
    R = Vasicek(t, mT, r0, k, rhat, sigma) # 100 results
    p = mean(exp(-R)) 
    return (1000 * p)
  }
  
  else if(style =="coupon"){
    period = (mT-t) * 2
    cf = rep(coupon, period)
    cf[period] = coupon + 1000
    R = array(NA, period)
    pv = array(NA, period)
    for(i in 1:period ){
      R[i] = mean(Vasicek(t = t, mT = i/2-t, r0 = r0, k = k, rhat = rhat, sigma = sigma)) # n rates
      pv[i] = cf[i] * exp( - R[i] )
    }
    return (sum(pv))
  }
}

A = function(t,mT){
  exp( (rhat - sigma^2/2/kappa^2) * (B(t,mT) - (mT-t)) - sigma^2/4/kappa*B(t,mT)^2 )
}

B = function(t,mT){
  1/kappa * (1 - exp(-kappa*(mT-t)))
}

# Explicit pure discount bond price P(t,T) with a par value of $1
bondprice.explicit = function(t,mT){
  rt = Vasicek.rt(t, mT, r0, k = kappa, rhat, sigma)
  A(t,mT) * exp(-B(t,mT) * rt) # returns a vector of prices
}

europeancall.1c = function(t, mT, S, strike, type = c("explicit","implicit")){
  if(type=="explicit"){
    # t<mT<S
    pTS = bondprice.explicit(t = mT, mT = S)
  }
  else if(type=="implicit"){
    pTS = bondprice(mT, S, r0, k=kappa, rhat, sigma, style = "zero")
  }
  R = Vasicek(t, mT, r0, k = kappa, rhat, sigma) # vector of discount rates
  call =  mean(exp(-R) * pmax(pTS * 1000 - strike, 0))
  return (call)
}

europeancall.1d = function(t, mT, S, strike){
  # P(T,S)
  period = S * 2
  cf = rep(30, period)
  cf[period] = 1030
  
  # r integration(T,S)
  # each time step is a day, assume 360 days a year
  dt = 1/360
  nstep = (S-mT) / dt
  r = matrix(NA, nrow = nSims, ncol = nstep) # 100 Monte Carlo simulations
  r[,1] = r0
  for( i in 2:nstep ){
    dw = rnorm(nSims)
    r[,i] = r[,i-1] + kappa * (rhat - r[,i-1]) * dt + sigma * sqrt(dt) * dw
  }
  integration = apply(r[,], 1, sum) * dt
  R = array(NA, period)
  pv = array(NA, period)
  for(i in 1:period ){
    R[i] = mean(Vasicek(mT, i/2-mT, r0 = r0, k = k, rhat = rhat, sigma = sigma)) # n rates
    pv[i] = cf[i] * exp( - R[i] )
  }
  pTS = sum(pv)
  
  R = Vasicek(t, mT, r0, k = kappa, rhat, sigma) # vector of discount rates
  call =  mean(exp(-R) * pmax(pTS - strike, 0))
  return (call)
}

#### COMPUTATION ####

# (a) zero-coupon bond
p_a = bondprice(t=0, mT = 0.5, r0 = r0, k = kappa, rhat = rhat, sigma = sigma, style="zero")

# (b) coupon paying bond
p_b = bondprice(t=0, mT = 4, r0 = r0, k = kappa, rhat = rhat, sigma = sigma, coupon = 30, style="coupon")

# (c) European call option on zero-coupon bond, monte carlo call pricing, explicit bond pricing
p_c = europeancall.1c(t=0, mT = 3/12, S = 0.5, strike = 980, type = "explicit") # S is NOT 3/12+0.5!!!!

# (d) European call on coupon paying bond. Monte Carlo for call and bond pricing.
p_d = europeancall.1d(t = 0, mT = 3/12, S = 4, strike = 980)

#############################
# Problem 2 CIR model, Monte Carlo, bond pricing
r0 = 0.05
sigma = 0.12
kappa = 0.92
rhat = 0.055

#### FUNCTION DEFINITION ####
CIR = function(t, mT, r0, k, rhat, sigma){
  # assume 360 days a year
  nstep = (mT-t) * 360
  # each time step is a day
  dt = 1
  r = matrix(NA, nrow = 100, ncol = nstep) # 100 Monte Carlo simulations
  r[,1] = r0
  for( i in 2:nstep ){
    dw = rnorm(100)
    r[,i] = r[,i-1] + kappa * (rhat - r[,i-1]) * dt + sigma * sqrt(r[,i-1]) * dw
  }
  integration = apply(r, 1, sum) # because dt = 1
  return( integration/360 ) # mind that rates are annual, convert to the rate over dt
}

bondprice.cir = function(t, mT, r0,k, rhat, sigma){
  R = CIR(t, mT, r0, k, rhat, sigma) # 100 results
  p = exp(-R)
  return (1000 * p)
}