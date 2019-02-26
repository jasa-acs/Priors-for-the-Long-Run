function r = bvarPLR(y,lags,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates the BVAR of Giannone, Lenza and Primiceri (2017)
%
% y:        data matrix
% lags:     number of lags in the VAR
% For additional options see setPLR.m
% Last modified: 07/01/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_of_args = nargin;


%% set BVAR priors (several options available, see setpriors.m)
%  if varargin=[] --> default settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setPLR;


%% data matrix manipulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimensions
[TT,n]=size(y);
k=n*lags+1;         % # coefficients for each equation


% constructing the matrix of regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=zeros(TT,k);
x(:,1)=1;
for i=1:lags
    x(:,1+(i-1)*n+1:1+i*n)=lag(y,i);
end

y0    = mean(y(1:lags,:),1);
x     = x(lags+1:end,:);
y     = y(lags+1:end,:);
[T,n] = size(y);


% MN prior mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=zeros(k,n);
diagb=ones(n,1);
diagb(pos)=0;   % Set to zero the prior mean on the first own lag for variables selected in the vector pos
b(2:n+1,:)=diag(diagb);

% starting values for the minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda0=.2;     % std of MN prior
Phi0=1;         % std of plr prior

% residual variance of AR(1) for each variable
SS=zeros(n,1);
for i=1:n
    ar1=ols1(y(2:end,i),[ones(T-1,1),y(1:end-1,i)]);
    SS(i)=ar1.sig2hatols;
end


if mn==1;
    inlambda=-log((MAX.lambda-lambda0)./(lambda0-MIN.lambda));
else
    inlambda=[];
end

nlr = sum(IndPhi>=1);

if SinglePhi==1
    nlr = min(nlr,1);
elseif SinglePhi ~=0
    error('SinglePhi should be either active (1) or OFF (0)')
end;

if size(HH,1) == n & length(IndPhi)<n
    error('When H is of full rank you must speficy the kind of prior for each column')
end;

if size(HH,1) < n &  length(IndPhi)<=size(HH,1)
   error('When H is or reduced rank you must specify not only the kind of prior for each colums of H but also for the orthogonal complement of H: flat (0) or invariant (1)')   
end;



if nlr>0   
    inPhi=-log((MAX.Phi-Phi0)/(Phi0-MIN.Phi))*ones(nlr,1);
else
    inPhi=[];
end



x0=[inlambda;inPhi];

%x0 = [inlambda; (1:length(inPhi))'];
if mn == 0 & sum(IndPhi==1)>0
    error('Cannot estimate the hyperparameters for the PLR since the marginal likelihood is equal to 0 if mn is not active')
end;


if length(x0)>0
    
H0=10*eye(length(x0));          % initial guess for the inverse Hessian


%% maximization of the posterior of the hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fh,xh,gh,H,itct,fcount,retcodeh] = csminwel('logMLbvarPLR',x0,H0,[],1e-10,1000,y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,HH,HHcompl,IndPhi,SinglePhi,y0,hyperpriors,priorcoef,HPmode);
r.postmax.itct=itct;                % #iteration before reaching maximum

else

xh = x0
    
end;

%% output of the maximization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAR coefficients and residuals at the posterior model
[fh,r.postmax.betahat,r.postmax.sigmahat,r.postmax.lambda,r.postmax.Phi,r.HH]=logMLbvarPLR(xh,y,x,lags,T,n,b,MIN,MAX,SS,Vc,pos,mn,HH,HHcompl,IndPhi,SinglePhi,y0,hyperpriors,priorcoef,HPmode);


r.lags = lags;                      % # lags
r.postmax.SSar1=SS;                 % residual variance of AR(1) for each variable
r.postmax.logPost=-fh;              % value of the posterior of the hyperparameters at the peak


    


%% forecasts at the posterior mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Fcast==1
    Y=[y;zeros(hz(end),n)];
    for tau=1:max(hz)
        xT=[1;reshape(Y([T+tau-1:-1:T+tau-lags],:)',k-1,1)]';
        Y(T+tau,:)=xT*r.postmax.betahat;
    end
    r.postmax.forecast=Y(T+hz,:);
end


r.IndPhi = IndPhi;