%% MGMT237G Project 6 Wei Wei
%% ********************************************
% Problem 1 Fixed Strike Lookback Call and Put
mT=1;
r=0.03;
S0=98;
X=100;
sigma=0.12:0.04:0.48;
step=1000;
path=1000;

% computation
dt=mT/step;
Call=zeros(numel(sigma),1);
Put=zeros(numel(sigma),1);
for i=1:numel(sigma)
    R = exp((r-(sigma(i)^2)/2)*dt+sigma(i)*sqrt(dt)*randn(path,step));
    St =horzcat(repmat(S0,path,1),S0*cumprod(R,2));
    Max=max(St,[],2);
    Min=min(St,[],2);
    Call(i)=mean(exp(-r*mT)*max(Max-X,0));
    Put(i)=mean(exp(-r*mT)*max(X-Min,0));
end

% graph output
figure
plot(sigma,Call);
xlabel('Volatility');
ylabel('Lookback Call Price');
title('Proj6_1a');

figure
plot(sigma,Put);
xlabel('Volatility');
ylabel('Lookback Put Price');
title('Proj6_1b');

%% **********************************************
% Problem 2 Processes with Jumps, loan default pricing
%% 2(i)
lambda1=0.2;
lambda2=0.4;
T=5;
[D, Prob, Et] = p2func(lambda1, lambda2, T)
%% 2(ii)
lambda1=0.2;
lambda2=0:0.1:0.8;
T=3:1:8;
D_a=zeros(numel(T),numel(lambda2));
Prob_a=zeros(numel(T),numel(lambda2));
Et_a=zeros(numel(T),numel(lambda2));
for i=1:numel(T)
    for j=1:numel(lambda2)
        [D_a(i,j), Prob_a(i,j), Et_a(i,j)] = p2func(lambda1, lambda2(j), T(i))
    end
end
x = zeros(numel(T),numel(lambda2));
for col=1:numel(lambda2)
    x(:,col)=T;
end
y = zeros(numel(T),numel(lambda2));
for row=1:numel(T)
    y(row,:)=lambda2;
end
figure
subplot(2,2,1);
surf(x,y,D_a);
xlabel('T=3:8');
ylabel('lambda2=0:0.8');
zlabel('Default Option Value');
title('option value');
subplot(2,2,2);
surf(x,y,Prob_a);
xlabel('T=3:8');
ylabel('lambda2=0:0.8');
zlabel('default probability');
title('default probability');
subplot(2,2,3);
surf(x,y,Et_a);
xlabel('T=3:8');
ylabel('lambda2=0:0.8');
zlabel('Et');
title('Expected Exercise Time');


lambda1=0.05:0.05:0.4;
lambda2=0.4;
D_b=zeros(numel(T),numel(lambda1));
Prob_b=zeros(numel(T),numel(lambda1));
Et_b=zeros(numel(T),numel(lambda1));
for i=1:numel(T)
    for j=1:numel(lambda1)
        [D_b(i,j), Prob_b(i,j), Et_b(i,j)] = p2func(lambda1(j), lambda2, T(i))
    end
end
x = zeros(numel(T),numel(lambda1));
for col=1:numel(lambda1)
    x(:,col)=T;
end
y = zeros(numel(T),numel(lambda1));
for row=1:numel(T)
    y(row,:)=lambda1;
end
figure
subplot(2,2,1);
surf(x,y,D_b);
xlabel('T=3:8');
ylabel('lambda1=0.05:0.4');
zlabel('Default Option Value');
title('option value');
subplot(2,2,2);
surf(x,y,Prob_b);
xlabel('T=3:8');
ylabel('lambda1=0.05:0.4');
zlabel('default probability');
title('default probability');
subplot(2,2,3);
surf(x,y,Et_b);
xlabel('T=3:8');
ylabel('lambda1=0.05:0.4');
zlabel('Et');
title('Expected Exercise Time');
