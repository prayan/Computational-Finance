function [D, Prob, Et]  =  p2func(lambda1, lambda2, mT)
    V0 = 20000;
    L0 = 22000;
    mu = -0.1;
    sigma = 0.2;
    gamma = -0.4;
    r0 = 0.02;
    delta = 0.25;
    alpha = 0.7;
    ypsilon = 0.95;
    step = 500;
    dt = mT/step;
    R = r0+delta*lambda2;
    r = R/12;
    n = 12*mT;
    PMT = L0*r/(1-1/(1+r).^n);
    a = PMT/r;
    b = PMT/r/(1+r).^n;
    c = 1+r;
    beta = (ypsilon-alpha)/mT;
    path = 5000;
    Vt = zeros(path,step+1);
    Vt(:,1) = V0;

    %% normal random variables
    Norm = randn(path,step);
    Lt = zeros(1,step+1);
    Lt(1) = L0;
    Lt(2:end) = a-b*c.^(12*dt*(1:step));
    Qt = alpha+beta*dt*(0:1:step);
    QLt = Qt.*Lt;
    
    %% Poisson Distribution
    U1 = rand(path,step);
    X1 = -log(U1)/lambda1;
    clear U1;
    Jt = cumsum(X1,2);
    clear X1;
    Jt = round(Jt*100)/100;
    U2 = rand(path,step);
    X2 = -log(U2)/lambda2;
    clear U2;
    Nt = cumsum(X2,2);%For calculating S
    clear X2;
    Nt = round(Nt*100)/100;
    S = Nt(:,1);
    
    %% Vt of different paths
    Index = zeros(path,step+1);
    Q = zeros(path,1);
    for j = 1:path
        temp = Jt(j,:);
        temp = temp(temp<= mT);
        for i = 2:step+1;
            t = (i-1)*dt;
            Vt(j,i) = Vt(j,i-1)+mu*Vt(j,i-1)*dt+sigma*Vt(j,i-1)*sqrt(dt).*...
            Norm(j,i-1)+gamma*(isscalar(find(temp == t,1,'first')));
            Index(j,i) = (Vt(j,i)<= QLt(i));
        end
        if isscalar(find(Index(j,:),1,'first'))
            Q(j) = (find(Index(j,:),1,'first')-1)*dt;
        else Q(j) = mT+1;% Assigning a larger-than-T value to Q that won't be used
        end
    end
    Tao = min(Q,S);
    
    Payoff = zeros(path,1);
    Default = 0;
    for j = 1:path 
        if and(Tao(j) == S(j),S(j)<= mT) == 1
            Locate = round(S(j)/dt+1);
            Payoff(j) = abs(Lt(Locate)-ypsilon*Vt(j,Locate));
        elseif and(Tao(j) == Q(j),Q(j)<= mT) == 1
            Locate = round(Q(j)/dt+1);
            Payoff(j) = max(Lt(Locate)-ypsilon*Vt(j,Locate),0);
        else Payoff(j) = 0; 
            Default = Default+1;
        end
    end
    D = mean(Payoff);
    Prob = (path-Default)/path;
    Et = mean(Tao(Tao<= mT));
end


