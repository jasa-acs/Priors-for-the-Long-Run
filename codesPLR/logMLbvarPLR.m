function [logML,betahat,sigmahat,lambda,Phi,HH]=logMLbvarPLR(par,y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,HH,HHcompl,IndPhi,SinglePhi,y0,hyperpriors,priorcoef,HPmode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the log-posterior (or the logML if hyperpriors=0),
% the posterior mode of the coefficients and the covariance matrix of the residuals of the BVAR of
% Giannone, Lenza and Primiceri (2017)
%
% Last modified: 10/13/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d=n+2;
cont = 1;
if mn == 1
    lambda=MIN.lambda+(MAX.lambda-MIN.lambda)/(1+exp(-par(cont)));
    cont = cont +1;
elseif mn == -1 %Dogmatic
    lambda = 0;
elseif mn == -2 %Mode
    lambda = HPmode.lambda;
elseif mn == 0  %OFF
    lambda = inf;
else
    error('mn can only take the values 0, 1, -1, -2')
end;

alpha=2;
psi=SS;
for jlr = 1:length(IndPhi)
    if IndPhi(jlr)==0 %OFF
        Phi(jlr) = inf;
    elseif IndPhi(jlr)==-1 %Dogmatic
        Phi(jlr) = 0.0001;
    elseif IndPhi(jlr)==-2 %Mode
        Phi(jlr) = HPmode.Phi;
    elseif IndPhi(jlr)==1 %Optimize
        Phi(jlr)=MIN.Phi+(MAX.Phi-MIN.Phi)./(1+exp(-par(cont)));
        if SinglePhi == 0
            cont = cont+1;
        end;
    else
        error('IndPhi can only take the values 0, 1, -1, -2')
    end    
end

%% setting up the priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1+n*lags;
omega=zeros(k,1);
omega(1)=Vc;

for i=1:lags
    omega(1+(i-1)*n+1:1+i*n)=(d-n-1)*(lambda^2)*(1/(i^alpha))./psi;
end

% prior scale matrix for the covariance of the shocks
PSI=diag(psi);

Td=0;
ydlr = [];
xdlr=[];
% dummy observations for the long run priors

rHH = size(HH,1); %rank of HH (HH must be full row rank)

if rHH<n
    HH = [HH;HHcompl];
end;

y0H = y0*HH';

HH_1 = inv(HH);

cont = 0;

for jlr = 1:length(IndPhi)
    if jlr <=rHH
        if IndPhi(jlr)~=0
            Phitemp = Phi(jlr);
            scale = y0H(jlr)/Phitemp;
            ydtemp = scale*HH_1(:,jlr)';
            xdtemp = [0  scale*repmat(HH_1(:,jlr)',1,lags)];
        end;
    elseif IndPhi(jlr)==1 | IndPhi(jlr)==-2
        Phitemp =Phi(jlr);
        JJ = rHH+1:n;
        ydtemp = (HH_1(:,JJ)*HH(JJ,:)*y0')'/Phitemp;
        xdtemp = [1/Phitemp repmat((HH_1(:,JJ)*HH(JJ,:)*y0')',1,lags)/Phitemp];
    elseif IndPhi(jlr)~=0
        error('IndPhi the null space can only take the values 0,1,-2')
    end;
    if IndPhi(jlr)~=0
        ydlr = [ydlr;ydtemp];
        xdlr = [xdlr;xdtemp];
    end;
    
end;


y=[y;ydlr];
x=[x;xdlr];
Td = Td+size(xdlr,1);




T=T+Td;

%% output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% posterior mode of the VAR coefficients
if mn == 0 %flat
    lambda = inf;
    betahat=(x'*x)\(x'*y); 
    epshat=y-x*betahat;    % VAR residuals
    sigmahat=(epshat'*epshat)/(T-k);% Posterior mode of the covariance matrix
    logML=-inf; % ML equal zero bacause of fla
end;

if mn~=0 
    if mn == -1 %Dogmatic
        betahat = b;
    else %informative
        betahat=(x'*x+diag(1./omega))\(x'*y+diag(1./omega)*b);
    end
    
    % VAR residuals
    epshat=y-x*betahat;
    
    % Posterior mode of the covariance matrix
    if mn ==-1
        sigmahat=(epshat'*epshat + PSI)/(T+d+n+1);
    else
        sigmahat=(epshat'*epshat + PSI + (betahat-b)'*diag(1./omega)*(betahat-b))/(T+d+n+1);
    end;
    
    % logML
    aaa=diag(sqrt(omega))*(x'*x)*diag(sqrt(omega));
    if mn ==-1
        bbb=diag(1./sqrt(psi))*(epshat'*epshat )*diag(1./sqrt(psi));
    else
        bbb=diag(1./sqrt(psi))*(epshat'*epshat + (betahat-b)'*diag(1./omega)*(betahat-b))*diag(1./sqrt(psi));
    end;
    
    eigaaa=real(eig(aaa)); eigaaa(eigaaa<1e-12)=0; eigaaa=eigaaa+1;
    eigbbb=real(eig(bbb)); eigbbb(eigbbb<1e-12)=0; eigbbb=eigbbb+1;
    
    logML = - n*T*log(pi)/2 + sum(gammaln((T+d-[0:n-1])/2)-gammaln((d-[0:n-1])/2)) +...
        - T*sum(log(psi))/2 - n*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
    
    
    if isempty(ydlr)~=1;
        yd=[ydlr];
        xd=[xdlr];
        
        % prior mode of the VAR coefficients
        % betahatd=(xd'*xd+diag(1./omega))\(xd'*yd+diag(1./omega)*b);
        betahatd=b;     % this is the case for our priors (the line above delivers the same but is numerically not very stable)
        
        % VAR residuals at the prior mode
        epshatd=yd-xd*betahatd;
            
    
        aaa=diag(sqrt(omega))*(xd'*xd)*diag(sqrt(omega));
        
        if mn ==-1
            bbb=diag(1./sqrt(psi))*(epshatd'*epshatd)*diag(1./sqrt(psi));
        else
            bbb=diag(1./sqrt(psi))*(epshatd'*epshatd + (betahatd-b)'*diag(1./omega)*(betahatd-b))*diag(1./sqrt(psi));
        end
        
        eigaaa=real(eig(aaa)); eigaaa(eigaaa<1e-12)=0; eigaaa=eigaaa+1;
        eigbbb=real(eig(bbb)); eigbbb(eigbbb<1e-12)=0; eigbbb=eigbbb+1;
        
        % normalizing constant
        norm = - n*Td*log(pi)/2 + sum(gammaln((Td+d-[0:n-1])/2)-gammaln((d-[0:n-1])/2)) +...
            - Td*sum(log(psi))/2 - n*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
        
        logML=logML-norm;
    end
    
    
    if hyperpriors==1;
        if mn == 1
            logML=logML+logPhipdf(lambda,priorcoef.lambda.k,priorcoef.lambda.theta);
        end
        
        temp = 0;
        for jlr = 1:length(IndPhi)
            if IndPhi(jlr) == 1 & temp==0
                logML=logML+logPhipdf(Phi(jlr),priorcoef.Phi.k,priorcoef.Phi.theta);
                if SinglePhi == 1
                    temp = 1;
                end
            end;
        end
    end;
    
    logML=-logML;
end;



function r=logPhipdf(x,k,theta);
r=(k-1)*log(x)-x/theta-k*log(theta)-gammaln(k);

function r=logIG2pdf(x,alpha,beta);
r=alpha*log(beta)-(alpha+1)*log(x)-beta./x-Philn(alpha);
