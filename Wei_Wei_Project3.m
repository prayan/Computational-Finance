% Wei Wei
% M237G Project 3 April 20 2016

%%%% Problem 1
dt = 0.001;
N = 10000;
x0 = 1;
y0 = 3/4;

Xt = zeros(N,2/dt+1);
Xt(:,1) = X0;
Yt = zeros(N,3/dt+1);
Yt(:,1) = Y0;

% X2 path
for i=1:2000
    t = dt * (i-1);
    dXt = (1/5 - 0.5 * Xt(:,1)) * dt + (2/3) * sqrt(t) * normrnd(0,1,N,1);
    Xt(:,i+1) = Xt(:,i) + dXt;
end
% Y3 path
for i=1:3000
    t = dt * (i-1);
    dYt = ( 2 / (1+t) * Yt(:,i) + (1 + t^3)/3 ) * dt + (1 + t^3) / 3 * sqrt(dt) * normrnd(0,1,N,1);
    Yt(:,i+1) = Yt(:,i) + dYt; 
end

i2 = 2001;
i3 = 3001;

% P(Y2>5)
Prob = sum(Yt(:,i2)>5)/N;
% E(X2^(1/3))
E1 = mean(sign(Xt(:,i2)).*abs(Xt(:,i2).^(1/3)));
% E(Y3)
E2 = mean(Yt(:,i3));
% E(X2Y2*1(PX2>1))
Z = Xt(:,i2)>1;
E3 = mean(Xt(:,i2).*Yt(:,i2).*Z);


%%%% Problem 2
dt = 0.001;
N = 10000;
X0 = 1;
Xt = zeros(N,3/dt+1);

for i=1:3000
    t = dt * (i-1);
    dXt = (1/4) * Xt(:,i) * dt + (1/3) * Xt(:,i) * sqrt(t).* normrnd(0,1,N,1) - (3/4) * Xt(:,i) * sqrt(t).* normrnd(0,1,N,1);
    Xt(:,i+1) = Xt(:,i) + dXt; 
end 

 Y3 = exp(-0.08*3 + 1/3 * sqrt(3) .* normrnd(0,1,N,1) + 3/4 * sqrt(3).* normrnd(0,1,N,1));
 E1 = mean(sign(1+Xt(:,i3)) .* abs((1+Xt(:,i3)) .^ (1/3)));
 E2 = mean((1+Y3) .^ (1/3));
 
 %%%% Problem 3
 %%%% 3(a)
S0=88;
r=0.04;
sigma=0.2;
T=5;
K=100;

N=10000;
F=normrnd(0,1,N,1);
% variance reduction
callplus = mean(max(S0*exp(sigma*sqrt(T)*F+(r-sigma^2/2)*T)-K,0)*exp(-r*T));
callminus = mean(max(S0*exp(sigma*sqrt(T)*(-F)+(r-sigma^2/2)*T)-K,0)*exp(-r*T));
C1=(callplus + callminus)/2;

%%%% 3(b)
S = 20;
T = 0.5;
K = 20;
r = 0.04;
sigma = 0.25;
%Computes European call option value in the Black-Scholes model
d1 = (log(S/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = (log(S/K)+(r-sigma^2/2)*T)/(sigma*sqrt(T));
C2 = S*normcdf(d1)-K*exp(-r*T)*normcdf(d2);

%%%% 3(c)
t = 0.5;
K = 20;
r = 0.04;
sigma = 0.25;
S0 = linspace(15,25,11);
D = zeros(1,11);
T = zeros(1,11);
V = zeros(1,11);
G = zeros(1,11);
R = zeros(1,11);
for i=1:11
    d1 = (log(S0(i)/K) + (r+1/2*sigma^2)*t) / (sigma*sqrt(t));
    d2 = (log(S0(i)/K) + (r-1/2*sigma^2)*t) / (sigma*sqrt(t));
    D(i) = normcdf(d1);
    T(i) = (-S0(i) * sigma * normpdf(d1)) / (2*sqrt(t)) - r*K*exp(-r*t) * normcdf(d2);
    V(i) = S0(i) * sqrt(t) * normpdf(d1);
    G(i) = 1 / (S0(i) * sigma * sqrt(t)) * normpdf(d1);
    R(i) = K * t * exp(-r*t) * normcdf(d2);
end
figure;
plot(S0,D,S0,T,S0,V,S0,G,S0,R);
legend('Delta','Theta','Vega','Gamma','Rho');
xlabel('S0');
ylabel('Greeks');


%%%% Problem 4
dt = 0.001;
N = 10000;
rho = -0.6;
r= 0.03;
S0 = 48;
V0 = 0.05;
sigma = 0.42;
alpha = 5.8;
beta = 0.0625;
K = 48;
T = 3;
% full truncation
StFull = zeros(N,T/dt+1);
VtFull = zeros(N,T/dt+1);
% partial truncation
StPartial = zeros(N,T/dt+1);
VtPartial = zeros(N,T/dt+1);
% reflection
StRefl = zeros(N,T/dt+1);
VtRefl = zeros(N,T/dt+1);
% Build matrixs for the maths of these three methods;
StFull(:,1) = S0;
VtFull(:,1) = V0;
StPartial(:,1) = S0;
VtPartial(:,1) = V0;
StRefl(:,1) = S0;
VtRefl(:,1) = V0;
Z1=normrnd(0,1,N,1);
Z2=normrnd(0,1,N,1);
for i=1:3000
    dWt1 = sqrt(dt)*Z1;
    dWt2 = sqrt(dt)*rho*Z1+sqrt(dt)*sqrt(1-rho^2)*Z2;
    dStF = r*StFull(:,i)*dt+sqrt(max(VtFull(:,i),0)).*StFull(:,i).*dWt1;
    dVtF = alpha*(beta-max(VtFull(:,i),0))*dt+sigma*sqrt(max(VtFull(:,i),0)).*dWt2;
    dStP = r*StPartial(:,i)*dt+sqrt(max(VtPartial(:,i),0)).*StPartial(:,i).*dWt1;
    dVtP = alpha*(beta-VtPartial(:,i))*dt+sigma*sqrt(max(VtPartial(:,i),0)).*dWt2;
    dStR = r*StRefl(:,i)*dt+sqrt(abs(VtRefl(:,i))).*StRefl(:,i).*dWt1;
    dVtR = alpha*(beta-abs(VtRefl(:,i)))*dt+sigma*sqrt(abs(VtRefl(:,i))).*dWt2;
    StFull(:,i+1) = StFull(:,i) + dStF;
    VtFull(:,i+1) = VtFull(:,i) + dVtF;
    StPartial(:,i+1) = StPartial(:,i) + dStP;
    VtPartial(:,i+1) = VtPartial(:,i) + dVtP;
    StRefl(:,i+1) = abs(StRefl(:,i)) + dStR;
    VtRefl(:,i+1) = abs(VtRefl(:,i)) + dVtR;
end
T_end = T/dt+1;
C1 = exp(-r*T) * mean(max(StFull(:,T_end) - K,0));
C2 = exp(-r*T) * mean(max(StPartial(:,T_end) - K,0));
C3 = exp(-r*T) * mean(max(StRefl(:,T_end) - K,0));


%%%% Problem 5
N = 100;
% 5(a)
U1 = unifrnd(0,1,N,2);
% 5(b)
U2 = zeros(N,2);
U2(:,1) = HaltonSeq(N,2);
U2(:,2) = HaltonSeq(N,7);
% 5(c)
U3=zeros(N,2);
U3(:,1) = HaltonSeq(N,2);
U3(:,2) = HaltonSeq(N,4);
% 5(d)
figure;
scatter(U1(:,1),U1(:,2));
figure;
scatter(U2(:,1),U2(:,2));
figure;
scatter(U3(:,1),U3(:,2));
% 5(e)
N = 10000;
% base (2,4)
H1 = zeros(N,2);
H1(:,1) = HaltonSeq(N,2);
H1(:,2) = HaltonSeq(N,4);
% base(2,7)
H2 = zeros(N,2);
H2(:,1) = HaltonSeq(N,2);
H2(:,2) = HaltonSeq(N,7);
% base(5,7)
H3 = zeros(N,3);
H3(:,1) = HaltonSeq(N,5);
H3(:,2) = HaltonSeq(N,7);

seq1 = sign(cos(2*pi*H1(:,2))) .* abs((cos(2*pi*H1(:,2))).^ (1/3));
Integral1 = mean(exp(-H1(:,1) .* H1(:,2)) .* (sin(6*pi*H1(:,1)) + seq1));

seq2 = sign(cos(2*pi*H2(:,2))) .* abs((cos(2*pi*H2(:,2))).^ (1/3));
Integral2 = mean(exp(-H2(:,1) .* H2(:,2)) .* (sin(6*pi*H2(:,1)) + seq2));

seq3 = sign(cos(2*pi*H3(:,2))) .* abs((cos(2*pi * H3(:,2))).^ (1/3));
Integral3 = mean(exp(-H3(:,1).* H3(:,2)) .* (sin(6*pi*H3(:,1)) + seq3));