# MGMT237G Project 7 Finite Differential Method
# Wei Wei

###############################
# Problem 1
###############################
currPrice = 10
k = 10 # strike price
mt = 0.5 # maturity
sigma = 0.2
rf = 0.04
dt = 0.002
dx = c( sigma*sqrt(dt), sigma*sqrt(3*dt), sigma*sqrt(4*dt) )
type = "put"
style = "european"

#######################################
# function definition

# (a) Explicit Finite-Difference Method
fde.log = function(s, k, r, t, sd,
                    n = ceiling(1000*t), m = 2*ceiling(sqrt(3*n)),
                    type = c("call", "put"), style = c("european", "american"),
                    grid = FALSE) {
  type = match.arg(type); style = match.arg(style)
  
  dt = t / n
  m = m + m%%2                         # Ensure m is even.
  ## Set stock price limits to +/- 3 standard deviations.
  z.lim = log(s) + 3*sd*sqrt(t)*c(min=-1, max=1)
  dz = unname(diff(z.lim)) / m
  z.seq = z.lim['min'] + 0:m*dz        # vector, m+1 elements
  
  f = matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m = function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call=pmax(exp(z.seq)-k,0), put=pmax(k-exp(z.seq),0))
  
  a = (1 + r*dt)^-1 * (-dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  b = (1 + r*dt)^-1 * (1 - dt/dz^2*sd^2)
  c = (1 + r*dt)^-1 * (dt/(2*dz)*(r - 1/2*sd^2) + dt/(2*dz^2)*sd^2)
  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    j.seq = rep(g2m(0:2), times=m-1) + rep(0:(m-2), each=3)
    f[i,g2m(1:(m-1))] = matrix(f[i+1,j.seq], ncol=3, byrow=TRUE) %*% c(a,b,c)
    
    if (type == 'call') {               # m: ∂C/∂S ≈ 1
      f[i,g2m(m)] = f[i,g2m(m-1)] + exp(z.seq[g2m(m)]) - exp(z.seq[g2m(m-1)])
      f[i,g2m(0)] = f[i,g2m(1)]        # 0: ∂C/∂S ≈ 0
    }
    else if (type == 'put') {           # m: ∂C/∂S ≈ 0
      f[i,g2m(m)] = f[i,g2m(m-1)]      # 0: ∂C/∂S ≈ 1
      f[i,g2m(0)] = f[i,g2m(1)] - (exp(z.seq[g2m(m)]) - exp(z.seq[g2m(m-1)]))
      if (style == 'american')
        f[i,] = pmax(f[i,], k - exp(z.seq))
    }
  }
  
  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

# (b) Implicit Finite-Difference Method
fdi.log = function(s, k, r, t, sd,
                    n = ceiling(1e3*t), m = 2*ceiling(sqrt(3*n)),
                    type = c("call", "put"), style = c("european", "american"),
                    grid = FALSE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  type = match.arg(type); style = match.arg(style)
  
  dt = t / n
  m = m + m%%2                         # Ensure m is even.
  ## Set stock price limits to +/- 3 standard deviations.
  z.lim = log(s) + 3*sd*sqrt(t)*c(min=-1, max=1)
  dz = unname(diff(z.lim)) / m
  z.seq = z.lim['min'] + 0:m*dz        # vector, m+1 elements
  
  f = matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m = function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call=pmax(exp(z.seq)-k,0), put=pmax(k-exp(z.seq),0))
  f[,g2m(m)] = switch(type, call=exp(z.seq[g2m(m)])-k, put=0)
  f[,g2m(0)] = switch(type, call=0,                    put=k-exp(z.seq[g2m(0)]))
  
  a = dt/(2*dz)*(r - 1/2*sd^2) - dt/(2*dz^2)*sd^2
  b = 1 + dt/dz^2*sd^2 + r*dt
  c = -dt/(2*dz)*(r - 1/2*sd^2) - dt/(2*dz^2)*sd^2
  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    A = matrix(0, nrow=m-1, ncol=m-1)
    B = f[i+1,g2m(1:(m-1))]
    for (j in 1:(m-1)) {
      if (j <= 1) {                     # j == 1
        A[j,1:2] = c(b, c)
        B[1] = B[1] - a*f[i,g2m(0)]
      }
      else if (j < m-1)
        A[j,(j-1):(j+1)] = c(a, b, c)
      else {                            # j == m-1
        A[j,(m-2):(m-1)] = c(a, b)
        B[m-1] = B[m-1] - c*f[i,g2m(m)]
      }
    }
    f[i,g2m(1:(m-1))] = solve(A, B)
    
    if (type == 'put' && style == 'american')
      f[i,] = pmax(f[i,], k - exp(z.seq))
  }
  
  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

# (c) Crank-Nicolson Finite-Difference method
fdcn.log = function(s, k, r, t, sd,
                     n = ceiling(1e3*t), m = 2*ceiling(sqrt(3*n)),
                     type = c("call", "put"), style = c("european", "american"),
                     grid = FALSE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  type = match.arg(type); style = match.arg(style)
  
  dt = t / n
  m = m + m%%2                         # Ensure m is even.
  ## Set stock price limits to +/- 3 standard deviations.
  z.lim = log(s) + 3*sd*sqrt(t)*c(min=-1, max=1)
  dz = unname(diff(z.lim)) / m
  z.seq = z.lim['min'] + 0:m*dz        # vector, m+1 elements
  
  f = matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m = function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call=pmax(exp(z.seq)-k,0), put=pmax(k-exp(z.seq),0))
  
  p.u = -dt/4*(sd^2/dz^2 + (r - 1/2*sd^2)/dz)
  p.m = 1 + dt*(sd^2/(2*dz^2) + r/2)
  p.d = -dt/4*(sd^2/dz^2 - (r - 1/2*sd^2)/dz)
  tridiag = rbind(0, cbind(diag(c(rep(p.u, m-1), 1)), 0)) +
    diag(c(1, rep(p.m, m-1), -1)) +
    rbind(cbind(0, diag(c(-1, rep(p.d, m-1)))), 0)
  c = tridiag
  d = diag(2, m+1) - tridiag           # Note: 1st & last rows corrected later.
  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    rhs = d %*% f[i+1,g2m(m:0)]
    if (type == 'call') {
      rhs[g2m(0)] = exp(z.seq[g2m(m)]) - exp(z.seq[g2m(m-1)])
      rhs[g2m(m)] = 0
    }
    else if (type == 'put') {
      rhs[g2m(0)] = 0
      rhs[g2m(m)] = exp(z.seq[g2m(m)]) - exp(z.seq[g2m(m-1)])
    }
    f[i,g2m(m:0)] = solve(c, rhs)
    
    if (type == 'put' && style == 'american')
      f[i,] = pmax(f[i,], k - exp(z.seq))
  }
  
  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

# Black-Scholes-Merton option pricing
bsm = function(s, k, r, t, sigma, type = c("call", "put")){
  d1 = (log(s/k) + (r + 1/2*sigma^2)*t) / (sigma*sqrt(t))
  d2 = d1 - sigma*sqrt(t)
  price = switch(type, 
                 call = s*pnorm(d1) - k*exp(-r*t)*pnorm(d2), 
                 put = k*exp(-r*t)*pnorm(d2,lower.tail=F)-s*pnorm(d1,lower.tail=F))
  return(price)
}

###########################################
# Computation

# (a) Explicit
Pa = fde.log(currPrice, k, rf, mt, sigma, type=type, style=style)
# (b) Implicit
Pb = fdi.log(currPrice, k, rf, mt, sigma, type=type, style=style)
# (c) Crank-Nicolson
Pc = fdcn.log(currPrice, k, rf, mt, sigma, type=type, style=style)
# comparison with Black-Scholes-Merton formula value
P_bsm = bsm(currPrice, k, rf, mt, sigma, type=type)


currPrice = seq(4,16,1)
P_errors.euro = matrix(NA, nrow = length(currPrice), ncol = 4)
P_errors.euro[,1] = currPrice
colnames(P_errors.euro) = c("currPrice", "explicit", "implicit", "Crank-Nicolson")
for(i in 1:length(currPrice)){
  pa = fde.log(currPrice[i], k, rf, mt, sigma, type=type, style=style)
  pb = fdi.log(currPrice[i], k, rf, mt, sigma, type=type, style=style)
  pc = fdcn.log(currPrice[i], k, rf, mt, sigma, type=type, style=style)
  pbsm = bsm(currPrice[i], k, rf, mt, sigma, type=type)
  P_errors.euro[i,2] = pa/pbsm - 1
  P_errors.euro[i,3] = pb/pbsm - 1
  P_errors.euro[i,4] = pc/pbsm - 1
}

###############################
# Problem 2
###############################

#######################################
# function definition

# (a) Explicit Finite-Difference method
fde = function(s, k, r, t, sd, n = ceiling(1e3*t), m = 2*ceiling(sqrt(3*n)),
                type = c("call", "put"), style = c("european", "american"),
                grid = FALSE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  type = match.arg(type); style = match.arg(style)
  
  dt = t / n
  m = m + m%%2                         # Ensure m is even.
  s.lim = c(max=2*s, min=0)
  ds = unname(s.lim['max'] - s.lim['min']) / m
  s.seq = s.lim['min'] + 0:m*ds        # vector, m+1 elements
  
  f = matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m = function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call = pmax(s.seq - k, 0), put = pmax(k - s.seq, 0))
  
  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    for (j in (m-1):1) {
      a = (1 + r*dt)^-1 * (-1/2*r*j*dt + 1/2*sd^2*j^2*dt)
      b = (1 + r*dt)^-1 * (1 - sd^2*j^2*dt)
      c = (1 + r*dt)^-1 * (1/2*r*j*dt + 1/2*sd^2*j^2*dt)
      f[i,g2m(j)] = t(c(a,b,c)) %*% f[i+1,g2m((j-1):(j+1))]
    }
    
    if (type == 'call') {               # m: ∂C/∂S ≈ 1
      f[i,g2m(m)] = f[i,g2m(m-1)] + s.seq[g2m(m)] - s.seq[g2m(m-1)]
      f[i,g2m(0)] = f[i,g2m(1)]        # 0: ∂C/∂S ≈ 0
    }
    else if (type == 'put') {           # m: ∂C/∂S ≈ 0
      f[i,g2m(m)] = f[i,g2m(m-1)]      # 0: ∂C/∂S ≈ 1
      f[i,g2m(0)] = f[i,g2m(1)] - (s.seq[g2m(1)] - s.seq[g2m(0)])
      if (style == 'american')
        f[i,] = pmax(f[i,], k - s.seq)
    }
  }
  
  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

# (b) implicit finite-difference method
fdi = function(s, k, r, t, sd, n = ceiling(1e3*t), m = 2*ceiling(sqrt(3*n)),
                type = c("call", "put"), style = c("european", "american"),
                grid = FALSE) {
  if (t <= 0) stop("t = ", t, " is nonpositive!")
  type = match.arg(type); style = match.arg(style)
  
  dt = t / n
  m = m + m%%2                         # Ensure m is even.
  s.lim = c(max=2*s, min=0)
  ds = unname(s.lim['max'] - s.lim['min']) / m
  s.seq = s.lim['min'] + 0:m*ds        # vector, m+1 elements
  
  f = matrix(rep(NA, (n+1)*(m+1)), nrow=n+1)
  g2m = function(i)  i + 1             # grid index to matrix index
  f[g2m(n),] = switch(type, call = pmax(s.seq - k, 0), put = pmax(k - s.seq, 0))
  f[,g2m(m)] = switch(type, call = s.seq[g2m(m)] - k,  put = 0)
  f[,g2m(0)] = switch(type, call = 0,                  put = k)
  
  for (i in g2m((n-1):0)) {             # Iterate from end to beginning.
    A = matrix(0, nrow=m-1, ncol=m-1)
    B = f[i+1,g2m(1:(m-1))]
    for (j in 1:(m-1)) {
      a = 1/2*r*j*dt - 1/2*sd^2*j^2*dt
      b = 1 + sd^2*j^2*dt + r*dt
      c = -1/2*r*j*dt - 1/2*sd^2*j^2*dt
      if (j <= 1) {                     # j == 1
        A[j,1:2] = c(b, c)
        B[1] = B[1] - a*f[i,g2m(0)]
      }
      else if (j < m-1)
        A[j,(j-1):(j+1)] = c(a, b, c)
      else {                            # j == m-1
        A[j,(m-2):(m-1)] = c(a, b)
        B[m-1] = B[m-1] - c*f[i,g2m(m)]
      }
    }
    f[i,g2m(1:(m-1))] = solve(A, B)
    
    if (type == 'put' && style == 'american')
      f[i,] = pmax(f[i,], k - s.seq)
  }
  
  if (grid) return(f) else return(f[g2m(0), g2m(m/2)])
}

#######################################
# Computation

ds = c(0.5, 1, 1.5)
type = "call"
style = "american"
currPrice = 10

# (a) explicit finite-difference American call
Ca = fde(currPrice, k, rf, mt, sigma, type=type, style=style)
# (b) implicit finite-difference American call
Cb = fdi(currPrice, k, rf, mt, sigma, type=type, style=style)
# (c) Crank-Nicolson American call
Cc = fdcn.log(currPrice, k, rf, mt, sigma, type=type, style=style)
# comparison with Black-Scholes-Merton formula value
C_bsm = bsm(currPrice, k, rf, mt, sigma, type=type)

type = "put"
Pa.amer = fde(currPrice, k, rf, mt, sigma, type=type, style=style)
Pb.amer = fdi(currPrice, k, rf, mt, sigma, type=type, style=style)
Pc.amer = fdcn.log(currPrice, k, rf, mt, sigma, type=type, style=style)
P_bsm.amer = bsm(currPrice, k, rf, mt, sigma, type=type)

# plot option price as function of current stock price
currPrice = seq(4,16,1)
plotAmericanOptions = function(type = c("call", "put")){
  prices = matrix(NA, nrow = length(currPrice), ncol = 3)
  for( i in 1:length(currPrice)){
    prices[i,1] = fde(currPrice[i], k, rf, mt, sigma, type=type, style=style)
    prices[i,2]= fdi(currPrice[i], k, rf, mt, sigma, type=type, style=style)
    prices[i,3] = fdcn.log(currPrice[i], k, rf, mt, sigma, type=type, style=style)
  }
  matplot(prices, type = c("l"), pch=1, col = 1:3, xlab = "current stock price", ylab = "option price")
  title( main = switch(type, call = "American call price as a function of s0",
                       put = "American put price as a function of s0"))
  legend("topleft", legend = c("currPrice", "explicit", "implicit", "Crank-Nicolson"), col=1:3, pch=1) # optional legend
}
plotAmericanOptions("call")
plotAmericanOptions("put")

