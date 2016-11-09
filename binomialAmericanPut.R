# binomial tree algorithm for American puts
# inputs should not be vectors
put.amer = function(s0, X, r, sigma, maturity){
  dt = 1/360
  step = maturity/dt
  R = exp(r*dt)
  u = 3/4
  d = 1/u
  p =  ( R-d ) / ( u-d ) # risk-neutral probability
  P = array(NA, step+1) # option price, each row is an observation
  for( i in 0:step ){
    P[i+1] = max( X - s0*u^(step-i)*d^i, 0)
  }
  for( j in (step-1):0 ){
    for( i in 0:j ){
      P[i+1] = max(( p*P[i+1] + (1-p) * P[i+2] )/R, X - s0*u^(j-i)*d^i)
    }
  }
  return (P[1])
}

r = 0.05
sigma = 0.3
X = 100
maturity = 1
s0 = seq(80,120,4)
puts = array(0, length(s0))
put = put.amer(s0[2],X,r,sigma,maturity)
for( i in 1:length(s0) ){
  puts[i] = put.amer(s0[i], X, r, sigma, maturity)
}
