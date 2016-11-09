% M237G Project 8 Fixed Income Wei Wei

%% Problem 1
%%%%%%%%%%%%%
clear;
clc;
r0 = 0.05;
sigma = 0.1;
kappa = 0.82;
rbar = 0.05;
%% 1(a) Vasicek model
%  1000 paths of r
dt = 1/365; % assume 365 days a year. each day is a time step
T = 0.5;
t = 0;
F = 1000;
paths = 1000;
N = round((T-t)/dt+1);
r_1 = zeros(paths,N);
r_1(:,1) = r0;
Norm = randn(paths,N-1);
for i = 2:N
    r_1(:,i) = r_1(:,(i-1))+kappa*(rbar-r_1(:,(i-1)))*dt+...
        sigma*sqrt(dt)*Norm(:,(i-1));
end
% bond pricing by Monte Carlo simulation
P_1a = mean(1000*exp(-sum(r_1,2)*dt))

%% 1£¨b£©
% 1000 paths of r
T = 4;
t = 0;
dt = 1/365;
N = round((T-t)/dt+1);
r_2 = zeros(paths,N);
r_2(:,1) = r0;
Norm = randn(paths,N-1);
for i = 2:N
    r_2(:,i) = r_2(:,(i-1))+kappa*(rbar-r_2(:,(i-1)))*dt+...
        sigma*sqrt(dt)*Norm(:,(i-1));
end
% coupon bond pricing
Time = floor((0.5:0.5:4)/dt+1 % Ti of cashflow
CF = [30*ones(1,7),1030];
for i = 1:numel(Time)
    R(:,i) = sum(r_2(:,1:Time(i)),2)*dt;
end
PVs = repmat(CF,paths,1).*exp(-R);
P_1b = mean(sum(PVs,2))

%% 1(c) explicit pure discount bond option pricing
S = 0.5;
T = 0.25;
K = 980;
% 1000 paths of r
N = round((S-0)/dt+1);
r_3 = zeros(paths,N);
r_3(:,1) = r0;
Norm = randn(paths,N-1);
for i = 2:N
    r_3(:,i) = r_3(:,(i-1))+kappa*(rbar-r_3(:,(i-1)))*dt+...
        sigma*sqrt(dt)*Norm(:,(i-1));
end
% Bond Price P(T,S) by Explicit Formula 
B = 1/kappa*(1-exp(-kappa*(S-T)));
A = exp((rbar-sigma^2/2/kappa^2)*(B-(S-T))-sigma^2/4/kappa*B^2);
P_1c = A*exp(-B*r_3(:,round((T-0)/dt+1)))*1000;

Call_1c = mean(exp(-sum(r_3,2)*dt).*max(P_1c-K,0))

%% 1(d)
% 1000 paths of r
T = 0.25;
paths = 1000;
N = round(T/dt+1);
r_4 = zeros(paths,N);
r_4(:,1) = r0;
Norm = randn(paths,N-1);
for i = 2:N
r_4(:,i) = r_4(:,(i-1))+kappa*(rbar-r_4(:,(i-1)))*dt+...
sigma*sqrt(dt)*Norm(:,(i-1));
end

% Simulation of bond price for each path
S = 4;
N2 = round((S-T)/dt+1);
P_TS = zeros(paths,1);
Ti = 0.5:0.5:4;
for i = 1:paths
      r_4temp = zeros(paths,N2);
      Norm_temp = randn(paths,N2-1);
      r_4temp(:,1) = r_4(i,end);
      for j = 2:N2
          r_4temp(:,j) = r_4temp(:,(j-1))+kappa*(rbar-r_4temp(:,(j-1)))*dt+...
            sigma*sqrt(dt)*Norm_temp(:,(j-1));
      end
      temp_S = zeros(8,1);
      for k = 1:8
         temp_S(k) = round((Ti(k)-T)/dt+1);
         R = sum(r_4temp(:,(1:temp_S(k))),2)*dt;
         P_TS(i) = P_TS(i)+CF(k)*mean(exp(-R));
      end
end

Call_1d = mean(exp(-sum(r_4,2)*dt).*max(P_TS-K,0))

%% 1(e)
% Solve for rstar
A_e = zeros(8,1);
B_e = zeros(8,1);
Ti = 0.5:0.5:4;
PI_e = 0;
syms x
for j = 1:8
    B_e(j) = 1/kappa*(1-exp(-kappa*(Ti(j)-T)));
    A_e(j) = exp((rbar-sigma^2/2/kappa^2)*(B_e(j)-(Ti(j)-T))-...
        sigma^2/4/kappa*B_e(j)^2);
    PI_e = PI_e+(A_e(j)*exp(-B_e(j)*x))*CF(j);
end
eqn =  PI_e =  = K;
rstar = solve(eqn,x);

% Calculating P(t,T) using explicit formula
Bt_T = 1/kappa*(1-exp(-kappa*(T-t)));
At_T = exp((rbar-sigma^2/2/kappa^2)*(Bt_T-(T-t))-sigma^2/4/kappa*Bt_T^2);
Pt_T = At_T*exp(-Bt_T*r0);
% Assgin Ki P(t,Ti) and find call prices
c_e = zeros(8,1);
Ki = zeros(8,1);
Pt_Ti = zeros(8,1);
for j = 1:8
    Ki(j) = A_e(j)*exp(-B_e(j)*rstar);
    B = 1/kappa*(1-exp(-kappa*(Ti(j)-t)));
    A = exp((rbar-sigma^2/2/kappa^2)*(B-(Ti(j)-t))-sigma^2/4/kappa*B^2);
    Pt_Ti(j) = A*exp(-B*r0);
    sigmap = sigma*(1-exp(-kappa*(Ti(j)-T)))/kappa*sqrt...
        ((1-exp(-2*kappa*(T-t)))/2/kappa);
    d1 = log(Pt_Ti(j)/Ki(j)/Pt_T)/sigmap+sigmap/2;
    d2 = d1-sigmap;
    c_e(j) = Pt_Ti(j).*normcdf(d1,0,1)-Ki(j)*Pt_T.*normcdf(d2,0,1);
end

Call_1e = CF*c_e

%% Problem 2
%%%%%%%%%%%%%%%
clear;
r0  =  0.05;
rbar  =  0.055;
kappa  =  0.92;
sigma  =  0.12;
T  =  .5;
S  =  1;
K  =  980;
dt  =  1/365;
L  =  1000;
t  =  0;
%% 2(a)
% Simulation of r from t to T
N  =  round((T-t)/dt+1);
paths  =  5000;
r_a  =  zeros(paths,N);
r_a(:,1)  =  r0;
Norm2a  =  randn(paths,N-1);
for i  =  2:N
    r_a(:,i)  =  r_a(:,i-1)+kappa*(rbar-r_a(:,i-1))*dt+sigma*sqrt...
        (r_a(:,i-1)*dt).*Norm2a(:,i-1);
end
R  =  sum(r_a(:,1:N),2)*dt;
% Simulation of P(T,S) using rT as initial interest rates
N2  =  round((S-T)/dt+1);
P_TS  =  zeros(paths,1);
for i  =  1: paths
  r_a2  =  zeros(paths/5,N2);
  r_a2(:,1)  =  r_a(i,end);
  Norm2a2  =  randn(paths/5,N2-1);
  for j  =  2:N2
      r_a2(:,j)  =  r_a2(:,j-1)+kappa*(rbar-r_a2(:,j-1))*dt+sigma*sqrt...
        (r_a2(:,j-1)).*Norm2a2(:,j-1)*sqrt(dt);
  end
  P_TS(i)  =  L*mean(exp(-sum(r_a2,2)*dt));
end
Call_2a  =  mean(exp(-R).*max(P_TS-K,0))

%% 2(b)
rmax  =  0.11;
rmin  =  0;
dr  =  0.0005;
dt  =  1/365;
r  =  (rmin:(dr):rmax)';
% Calculating P(T,S)
h1  =  sqrt(kappa^2+2*sigma^2);
h2  =  (kappa+h1)/2;
h3  =  2*kappa*rbar/sigma^2;
A_TS  =  (h1*exp(h2*(S-T))/(h2*(exp(h1*(S-T))-1)+h1))^h3;
B_TS  =  (exp(h1*(S-T))-1)/(h2*(exp(h1*(S-T))-1)+h1);
P_TS  =  L*A_TS*exp(-B_TS*r);
% Calculating option prices
F  =  zeros(numel(r),round(T/dt+1));
F(:,end)  =  max(P_TS-K,0);
PU  =  -0.5*dt*(kappa*(rbar-r)/dr+sigma^2*r/dr^2);
PM  =  1+dt*r+dt*sigma^2*r/dr^2;
PD  =  -0.5*dt*(-kappa*(rbar-r)/dr+sigma^2*r/dr^2);
A  =  zeros(numel(r));
A(1,1:2)  =  [1,-1];
A(end,(end-1:end))  =  [1,-1];
for j  =  2:(numel(r)-1)
       A(j,(j-1):(j+1))  =  [PU(j),PM(j),PD(j)];
end
temp  =  P_TS(1)-P_TS(2);
for i  =  round((T/dt)):-1:1
    B  =  F(:,i+1);
    B(1)  =  temp;
    B(end)  =  0;
    F(:,i)  =  A\B;
end
Call  =  horzcat(r,F(:,1));
Call_2b  =  Call(101,2)% Find the price when r  =  r0  =  0.05
%% 2(c)
% P(t,S)
A_tS  =  (h1*exp(h2*(S-t))/(h2*(exp(h1*(S-t))-1)+h1))^h3;
B_tS  =  (exp(h1*(S-t))-1)/(h2*(exp(h1*(S-t))-1)+h1);
P_tS  =  A_tS*exp(-B_tS*r0);
% P(t,T)
A_tT  =  (h1*exp(h2*(T-t))/(h2*(exp(h1*(T-t))-1)+h1))^h3;
B_tT  =  (exp(h1*(T-t))-1)/(h2*(exp(h1*(T-t))-1)+h1);
P_tT  =  A_tT*exp(-B_tT*r0);
% c(t,T,S)
theta  =  sqrt(kappa^2+2*sigma^2);
fai  =  2*theta/(sigma^2*exp(theta*(T-t)-1));
psi  =  (kappa+theta)/sigma^2;
A_TS  =  (h1*exp(h2*(S-T))/(h2*(exp(h1*(S-T))-1)+h1))^h3;
B_TS  =  (exp(h1*(S-T))-1)/(h2*(exp(h1*(S-T))-1)+h1);
rstar  =  log(L*A_TS/K)/B_TS;
Call_2c  =  L*P_tS*normcdf(2*rstar*(fai+psi+B_TS),4*kappa*rbar/sigma^2,2*fai^2*r0*exp(theta*(T-t))/(fai+psi+B_TS))-K*P_tT*normcdf(2*rstar*(fai+psi),4*kappa*rbar/sigma^2,2*fai^2*r0*exp(theta*(T-t))/(fai+psi))

%% Problem 3
%%%%%%%%%%%%%%%%%%
clear;
x0  =  0;
y0  =  0;
fai0  =  0.03;
r0  =  0.03;
rho  =  0.7;
a  =  0.1;
b  =  0.3;
sigma  =  0.03;
eta  =  0.08;
fait  =  0.03;
K  =  950;
L  =  1000;
T  =  0.5;
t  =  0;
S  =  1;
%% 3(a) Monte Carlo Simulation
paths  =  5000;
dt  =  1/365;
N  =  round(T/dt+1);
% Simulate two correlated normal series
N1  =  randn(paths,N-1);
N2  =  randn(paths,N-1);
Norm1  =  N1;
Norm2  =  N1*rho+N2*sqrt(1-rho^2);
clear N1 N2;
% Simulate r
xt  =  zeros(paths,N);
yt  =  zeros(paths,N);
r  =  zeros(paths,N);
r(:,1)  =  r0;
for i  =  2:N
    xt(:,i)  =  xt(:,i-1)-a*xt(:,i-1)*dt+sigma*sqrt(dt)*Norm1(:,i-1);
    yt(:,i)  =  yt(:,i-1)-b*yt(:,i-1)*dt+eta*sqrt(dt)*Norm2(:,i-1);
    r(:,i)  =  xt(:,i)+yt(:,i)+fait;
end
% Simulating P_TS Bond Prices for every path
P_TS  =  zeros(paths,1);
N2  =  round((S-T)/dt+1);
for i  =  1:paths
    xt_temp  =  zeros(paths/5,N2);
    xt_temp(:,1)  =  xt(i,end);
    yt_temp  =  zeros(paths/5,N2);
    yt_temp(:,1)  =  yt(i,end);
    rt_temp  =  zeros(paths/5,N2);
    rt_temp(:,1)  =  r(i,end);
    N3  =  randn(paths/5,N2-1);
    N4  =  randn(paths/5,N2-1);
    Norm3  =  N3;
    Norm4  =  N3*rho+N4*sqrt(1-rho^2);
    clear N3 N4
    for j  =  2:N2
    xt_temp(:,j)  =  xt_temp(:,j-1)-a*xt_temp(:,j-1)*dt+sigma*sqrt(dt)*...
        Norm3(:,j-1);
    yt_temp(:,j)  =  yt_temp(:,j-1)-b*yt_temp(:,j-1)*dt+eta*sqrt(dt)*...
        Norm4(:,j-1);
    rt_temp(:,j)  =  xt_temp(:,j)+yt_temp(:,j)+fait;
    end
    P_TS(i)  =  L*mean(exp(-sum(rt_temp,2)*dt));
end
        
% Put Prices for every path
R2  =  sum(r,2)*dt;
c  =  exp(-R2).*max(K-P_TS,0);
Put_3a  =  mean(c)
%% (b) Expicit Formula
% P(t,S)
V_tS  =  sigma^2/a^2*((S-t)+2/a*exp(-a*(S-t))-1/2/a*exp(-2*a*(S-t))-3/2/a)+...
    eta^2/b^2*((S-t)+2/b*exp(-b*(S-t))-1/2/b*exp(-2*b*(S-t))-3/2/b)+...
    2*rho*sigma*eta/a/b*((S-t)+(exp(-a*(S-t))-1)/a+(exp(-b*(S-t))-1)/b-...
    (exp(-(a+b)*(S-t))-1)/(a+b));
P_tS  =  exp(-S*fait-(1-exp(-a*(S-t)))/a*x0-(1-exp(-b*(S-t)))/b*y0+0.5*V_tS);
% P(t,T)
V_tT  =  sigma^2/a^2*((T-t)+2/a*exp(-a*(T-t))-1/2/a*exp(-2*a*(T-t))-3/2/a)+...
    eta^2/b^2*((T-t)+2/b*exp(-b*(T-t))-1/2/b*exp(-2*b*(T-t))-3/2/b)+...
    2*rho*sigma*eta/a/b*((T-t)+(exp(-a*(T-t))-1)/a+(exp(-b*(T-t))-1)/b-...
    (exp(-(a+b)*(T-t))-1)/(a+b));
P_tT  =  exp(-T*fait-(1-exp(-a*(T-t)))/a*x0-(1-exp(-b*(T-t)))/b*y0+0.5*V_tT);
% Delta
Delta_sq  =  sigma^2/2/a^3*(1-exp(-a*(S-T)))^2*(1-exp(-2*a*(T-t)))+...
    eta^2/2/b^3*(1-exp(-b*(S-T)))^2*(1-exp(-2*b*(T-t)))+...
    2*rho*sigma*eta/a/b/(a+b)*(1-exp(-a*(S-T)))*(1-exp(-b*(S-T)))*...
    (1-exp(-(a+b)*(T-t)));
Delta  =  sqrt(Delta_sq);
% Euro Put
Put_3b  =  -L*P_tS*normcdf(log(K*P_tT/L/P_tS)/Delta-0.5*Delta)+P_tT*K*...
    normcdf(log(K*P_tT/L/P_tS)/Delta+0.5*Delta)