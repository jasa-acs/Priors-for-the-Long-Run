% setPLR.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file sets up the default choices for the priors of the BVAR of 
% Giannone, Lenza and Primiceri (2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS:
%
% hyperpriors:   0 = no priors on hyperparameters
%                1 = reference priors on hyperparameters (default)
%
% Vc:            prior variance in the MN prior for the coefficients multiplying
%                the contant term (Default: Vc=10e6)
%           
% pos:           position of the variables that enter the VAR in first
%                differences and for which one might want to set the prior mean 
%                on the coefficient on the first own lag in the MN prior and the
%                prior mean of the sum-of-coefficients prior to 0 (instead of 1)
%                (Default: pos=[])
%
%
%
% fcast:         0 = does not generate forecasts at the posterior mode
%                1 = generates forecasts at the posterior mode (default)
%
% hz:            longest horizon at which the code generates forecasts
%                (default: maxhz=40)
%
% mn:            0 = Minnesota Prior is is OFF
%                1 = Estimate the overall tightness (default)
%               -1 = Dogmatically impose the prior
%               -2 = Sets the overall shrinkage at the mode of the hyperprior (.2)
%               
% HH:            Matrix to set the long run prior, each colums takes combinations of data
%                It should be a full row rank matrix
%
% HHcompl:       Completes HH is it not full rank
%                If not specified the codes construct it as the null space
%             
%               
%IndPhi:       Vector specifying the treatment of the PLR for each row of HH
%               0 = Prior is OFF 
%               1 = Estimate the LR tightness (default)
%              -1 = Dogmatically impose the prior
%              -2 = Sets the overall shrinkage at the mode of the hyperprior (1)
%               (if some raws of HH are not specified, the last index of IndPhi 
%               indicates the prior on HHcompl (is provided) or on the nulls space
%                        1: Invariant (default)
%                        0: flat 
%                       -2: Hyperprior Mode (2)
%                       
%SinglePhi:   1 = estimates the same degree of shrinkage for all active rows of HH
%             0 = Separate shrinkage for each active row (default)
%    
%Last modified: 10/13/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T,n] = size(y);
% Main options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if num_of_args>2
    for ii=1:2:size(varargin,2)
        eval([varargin{ii},'=',mat2str(varargin{ii+1})])
    end
end

%------------------------------
if ~exist('hyperpriors','var')
    hyperpriors=1;
end
r.setpriors.hyperpriors=hyperpriors;

%------------------------------
if ~exist('Vc','var')
    Vc=10e6;
end
r.setpriors.Vc=Vc;

%------------------------------
if ~exist('pos','var')
    pos=[];
end
r.setpriors.pos=pos;

%------------------------------
if ~exist('mn','var')
    mn=1;
end
r.setpriors.mn=mn;



%Priors for the long run
%------------------------------
if ~exist('HH','var')
    HH = [];
end
r.setpriors.HH=HH;


%Set of Hyperparameters
%------------------------------
if ~exist('IndPhi','var') | isempty(IndPhi)
    if isempty(HH)
        IndPhi = zeros(1,n);
        HH = eye(n);
    else
        temp = min(size(HH,1)+1,n);
        IndPhi = ones(1,temp);
    end;
end
r.setpriors.IndPhi=IndPhi;

%------------------------------
if ~exist('SinglePhi','var') | isempty(SinglePhi)
    SinglePhi=0;
end
r.setpriors.SinglePhi=SinglePhi;

%------------------------------
if ~exist('HHcompl','var') | isempty(HHcompl)
    HHcompl = null(HH)';
end



%------------------------------
if ~exist('Fcast','var')
    Fcast=1;
end
r.setpriors.Fcast=Fcast;

%------------------------------
if ~exist('hz','var')
    hz=[1:40]; 
else
    hz=[1:hz];
end
r.setpriors.hz=hz;



%% Other options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of the hyperpriors, if choosen

HPmode.lambda=.2;     % hyperpriors modes
HPmode.Phi=1;

HPsd.lambda=.4;       % hyperpriors std
HPsd.Phi = 1;


if hyperpriors==1;
    priorcoef.lambda=GammaCoef(HPmode.lambda,HPsd.lambda,0);  % coefficients of hyperpriors
    priorcoef.Phi=GammaCoef(HPmode.Phi,HPsd.Phi,0);
    
else
    priorcoef=[];
end

% bounds for maximization
MIN.lambda = 0.0001;
MIN.Phi    = 0.00001;
MAX.lambda = 5;
MAX.Phi  = 50;