# M237G Project 9 MBS, prepayment risk
# Wei Wei

nSims = 100
################
# Problem 1 ####
################

# FUNCTION DEFINITION #####

# Numerix Prepayment Model
NPM = function(mT, WAC, r0, kappa, rbar, sigma, L, OAS){
  N = 12 * mT # num of months, i.e. num of payments
  SY = c(0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22, 1.23, 0.98) # seasonality
  CPR = array(NA, N+1) # conditional prepayment rate
  TPP = array(NA, N) # total principal pmt
  IP = array(NA, N) # interest pmt
  Ct = array(NA, N) # cashflow at t
  PV = array(NA, N+1)
  PV[1] = L
  r = WAC / 12
  
  dt = 1/360
  rt = CIR.paths(r0, kappa, rbar, sigma, mT+10)
  for( i in 1:N ){
    rt10 = 1/10 * mean(apply(rt[,(i*30):(i*30+30)] * dt, 1, sum))
    RI = 0.28 + 0.14 * atan(-8.57 + 430 * (12*r - rt10))
    BU = 0.3 + 0.7 * PV[i] / L
    SG = min(1, i/30)
    SYt = SY[i%%12 + 12*(i%%12==0)]
    CPR[i+1] = RI * BU * SG * SYt
    TPP[i] = PV[i] * r * (1/(1-(1+r)^(-N+i-1))-1)+(PV[i]-PV[i]*r*(1/(1-(1+r)^(-N+i-1))-1))*(1-(1-CPR[i+1])^(1/12))
    IP[i] = PV[i] * r
    Ct[i] = TPP[i] + IP[i]
    PV[i+1] = PV[i] - TPP[i]
  }
  
  BondPrice = array(NA, N)
  h1 = sqrt(kappa^2 + 2*sigma^2)
  h2 = (kappa + h1)/2
  h3 = 2*kappa*rbar/sigma^2
  for( i in 1:N ){
    t = i/12 # years
    A = ( h1 * exp(h2*t) / (h2*(exp(h1*t)-1) + h1) )^h3
    B = ( exp(h1*t)-1 ) / (h2*(exp(h1*t)-1) + h1)
    BondPrice[i] = A * exp( -B * (r0+OAS) )
  }
  MBS = sum(Ct * BondPrice)
  return (MBS)
}

CIR.paths = function(r0, kappa, rbar, sigma, t){
  dt = 1/360 # dt = 1 day
  nStep = t / dt
  r = matrix(NA, nrow = nSims, ncol = nStep+1)
  r[,1] = r0
  for( i in 1:nStep ){
    dW = sqrt(dt) * rnorm(nSims)
    r[,i+1] = r[,i] + kappa * (rbar-r[,i]) * dt + sigma * sqrt(r[,i]) * dW
  }
  return (r)
}

# COMPUTATION #####

# 1(a)

"Input MBS parameters (separted by ,)"
input = readline(prompt = "Input MBS parameters: mT, WAC, r0, kappa, rbar, sigma, L, OAS: ")
input = as.numeric( unlist( strsplit(input, ",") ) )
p1.a = NPM(mT = input[1], WAC=input[2], r0=input[3], kappa=input[4], rbar=input[5], sigma=input[6], L=input[7], OAS=input[8])

# 1(b)
mT = 30
WAC = 0.08
r0 = 0.078
rbar = 0.08
sigma = 0.12
L = 100000
OAS = 0
kappa = seq(0.3,0.9,0.1)
p1.b = array(NA, dim = length(kappa))
for( i in 1:length(kappa) ){
  p1.b[i] = NPM(mT, WAC, r0, kappa[i], rbar, sigma, L, OAS)
}
plot(x = kappa, y = p1.b, type = "l", xlab = "kappa", ylab = "loan price")

# 1(c)
kappa = 0.6
rbar = seq(0.03,0.09,0.01)
p1.c = array(NA, dim = length(rbar))
for( i in 1:length(rbar) ){
  p1.c[i] = NPM(mT, WAC, r0, kappa, rbar[i], sigma, L, OAS)
}
plot(x = rbar, y = p1.c, type = "l", xlab = "rbar", ylab = "loan price")

################
# Problem 2 ####
################
# FUNCTION DEFINITION #####
# PSA Model of prepayments
PSA = function(mT, WAC, r0, kappa, rbar, sigma, L){
  N = 12 * mT # num of months, i.e. num of payments
  
  dt = 1/360
  rt = CIR.paths(r0, kappa, rbar, sigma, mT)
  BondPrice = array(NA, N)
  for( i in 1:N ){
    integration = apply(rt[,1:(i*30)]*dt, 1, sum)
    BondPrice[i] = mean( exp(-integration) )
  }
  
  CPR = array(0, N+1)
  TPP = array(0, N)
  IP = array(0, N)
  Ct = array(0, N)
  PV = array(0, N+1)
  PV[1] = L
  r = WAC/12
  
  for( i in 1:N ){
    if(i<=30){
      CPR[i+1] = CPR[i] + 0.002
    }
    else{
      CPR[i+1] = CPR[i]
    }
    TPP[i] = PV[i] * r * (1/(1-(1+r)^(-N+i-1)) - 1) +
      ( PV[i] - PV[i] * r *( 1 / ( 1 - (1+r)^(-N+i-1) ) - 1 ) ) * (1-(1-CPR[i+1])^(1/12))
    IP[i] = PV[i] * r
    Ct[i] = TPP[i] + IP[i]
    PV[i+1] = PV[i] - TPP[i]
  }
  
  Price = sum(Ct * BondPrice)
  return(Price)
}

# COMPUTATION #####

# 2(a)

"Input PSA MBS parameters (separted by ,)"
input2 = readline(prompt = "mT, WAC, r0, kappa, rbar, sigma, L: ")
input2 = as.numeric( unlist( strsplit(input2, ",") ) )
p2.a = PSA(mT = input2[1], WAC = input2[2], r0 = input2[3], kappa = input2[4], rbar = input2[5], sigma = input2[6], L = input2[7])

# 2(b)
kappa = seq(0.3, 0.9, 0.1)
mT = 30
WAC = 0.08
r0 = 0.078
rbar = 0.08
sigma = 0.12
L = 100000
p2.b = array(NA, length(kappa))
for( i in 1:length(kappa) ){
  p2.b[i] = PSA(mT, WAC, r0, kappa[i], rbar, sigma, L)
}
plot(x = kappa, y = p2.b, type = "l", xlab = "kappa", ylab = "loan price", main = "2(b)")

################
# Problem 3 Solve for OAS
################
Price = 110000 # market price of the loan
spread = function(OAS){
  abs(NPM(mT, WAC, r0, kappa, rbar, sigma, L, OAS) - Price)
}
test = optim(c(0.0030), fn = spread)
OAS = test$par # OAS = -0.05

################
# Problem 4 ####
################
x = OAS
y = 0.0005
p = NPM(mT, WAC, r0, kappa, rbar, sigma, L, x)
p.plus = NPM(mT, WAC, r0, kappa, rbar, sigma, L, x+y)
p.minus = NPM(mT, WAC, r0, kappa, rbar, sigma, L, x-y)
Duration = (p.minus - p.plus) / (2*y*p)
Convexity = (p.plus+p.minus-2*p) / (2*p*y^2)

################
# Problem 5 ####
################
# FUNCTION DEFINITION ####
NPM.IOPO = function(mT, WAC, r0, kappa, rbar, sigma, L, OAS){
  N = 12 * mT # num of months, i.e. num of payments
  SY = c(0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22, 1.23, 0.98) # seasonality
  CPR = array(NA, N+1) # conditional prepayment rate
  TPP = array(NA, N) # total principal pmt
  IP = array(NA, N) # interest pmt
  Ct = array(NA, N) # cashflow at t
  PV = array(NA, N+1)
  PV[1] = L
  r = WAC / 12
  
  dt = 1/360
  rt = CIR.paths(r0, kappa, rbar, sigma, mT+10)
  for( i in 1:N ){
    rt10 = 1/10 * mean(apply(rt[,(i*30):(i*30+30)] * dt, 1, sum))
    RI = 0.28 + 0.14 * atan(-8.57 + 430 * (12*r - rt10))
    BU = 0.3 + 0.7 * PV[i] / L
    SG = min(1, i/30)
    SYt = SY[i%%12 + 12*(i%%12==0)]
    CPR[i+1] = RI * BU * SG * SYt
    TPP[i] = PV[i] * r * (1/(1-(1+r)^(-N+i-1))-1)+(PV[i]-PV[i]*r*(1/(1-(1+r)^(-N+i-1))-1))*(1-(1-CPR[i+1])^(1/12))
    IP[i] = PV[i] * r
    Ct[i] = TPP[i] + IP[i]
    PV[i+1] = PV[i] - TPP[i]
  }
  
  BondPrice = array(NA, N)
  h1 = sqrt(kappa^2 + 2*sigma^2)
  h2 = (kappa + h1)/2
  h3 = 2*kappa*rbar/sigma^2
  for( i in 1:N ){
    t = i/12 # years
    A = ( h1 * exp(h2*t) / (h2*(exp(h1*t)-1) + h1) )^h3
    B = ( exp(h1*t)-1 ) / (h2*(exp(h1*t)-1) + h1)
    BondPrice[i] = A * exp( -B * (r0+OAS) )
  }
  IO = sum(IP * BondPrice)
  PO = sum(TPP * BondPrice)
  return (c(IO, PO))
}

# COMPUTATION ####
mT = 30
WAC = 0.08
r0 = 0.078
kappa = 0.6
sigma = 0.12
L = 100000
OAS = 0
rbar = seq(0.03, 0.09, 0.01)
IOPO = matrix(NA, ncol = length(rbar), nrow = 2)
for( i in 1:length(rbar) ){
  IOPO[,i] = NPM.IOPO(mT, WAC, r0, kappa, rbar[i], sigma, L, OAS)
}
IO = IOPO[1,]
PO = IOPO[2,]