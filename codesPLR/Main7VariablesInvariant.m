% Replcation files for the paper:
% "Priors for the Long-Run", by D. Giannone, G. Primiceri and M. Lenza
% 
% This code replicate the mean squared forecast errors in models with seven
% variables (figures 6.2)
clear;
clc;

% Loading the data
load y 
variables = {'Y','C','I','H','W','\pi','R'};
[T,n] = size(y);

addpath(['./subroutines']) %The folder with subroutines

lags = 5;  % n. of lags on the VAR
hz = 1:40; % Forecasting horizons in quarters
ini = 81;  % The beginning of the recursive estimation (out-of-sample)

%Inizializes the matrix to store the forecasts and the Mean Square Forecast Errrors
Yfcast=zeros(T,length(hz),n); 
MSFE=zeros(length(hz),n);    
MSFELC=zeros(length(hz),n);   

% Defines the matrix with the linear transformations used to elicit the Long Run Prior
%         Y   C   I   H   W  \pi  R
Ctr = [
          1   1   1   0   1   0   0;  %Y+C+I+W
         -1   1   0   0   0   0   0;  %C-Y
         -1   0   1   0   0   0   0;  %I-Y
          0   0   0   1   0   0   0;  %H
         -1   0   0   0   1   0   0;  %W-Y
          0   0   0   0   0   1   1;  %\pi + R
          0   0   0   0   0  -1   1]; %R-\pi 

transformations = {'Y+C+I+W','C-Y','I-Y','H','W-Y','R+\pi','R-\pi'};


% Define the different models
% There are several specifications you can run
% Just uncomment the text assoiciated to one of the specifications.
% You can run only one specification at a time and the line for the MSFEs will be added to the chart

% Different models are obtained by setting the following options (see setprior.m for details)
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
%               indicates the prior on the nulls spcae
%                        1: Invariant (default)
%                        0: flat 
%                       -2: Hyperprior Mode (2)
%                       
%SinglePhi:   1 = estimates the same degree of shrinkage for all active rows of HH
%             0 = Separate shrinkage for each active row (default)


%PLR benchmark
 ModLab = 'PLR';
 HH = Ctr;
 HHcompl = [];
 IndPhi = [1 1 1 1 1 1 1];
 SinglePhi = 0;
 mn = 1;
 color='r'; LS='-'; LW=3;

% %PLR Invariant
%  ModLab = 'PLR invariant';
%  HH = Ctr([1 6],:);
%  HHcompl = Ctr([2 3 4 5 7],:);
%  IndPhi = [1 1 1 ];
%  SinglePhi = 0;
%  mn = 1;
%  color= [1,.75,0]; LS=':'; LW=3;

% %PLR Invariant (except C-Y)
%  ModLab = 'PLR invariante (except C-Y)';
%  HH = Ctr([1 6 2],:);
%  HHcompl = Ctr([3 4 5 7],:);
%  IndPhi = [1 1 1 1];
%  SinglePhi = 0;
%  mn = 1;
%  color='r'; LS='-'; LW=3;
%  color= [.5,0,0]; LS='--'; LW=2;



% generating forecasts recursively
for t=ini:T 
    t
    pause(.5);
    r = bvarPLR(y(1:t,:),lags,'HH',HH,'mn',mn,'IndGamma',IndPhi,'SinglePhi',SinglePhi,'HHcompl',HHcompl);
    Yfcast(t,:,:)=r.postmax.forecast(:,:);
end

% computing the MSFE
for h=hz
        
    DYfcast=(squeeze(Yfcast(ini-h+max(hz):end-h,h,:))-y(ini-h+max(hz):end-h,:));%/h;
    DY=(y(ini+max(hz):end,:)-y(ini-h+max(hz):end-h,:));%/h;
    DYfcastDY=DYfcast-DY;
    
    MSFE(h,:)=mean((DYfcastDY).^2);
    
    DYfcastLC=(squeeze(Yfcast(ini-h+max(hz):end-h,h,:))-y(ini-h+max(hz):end-h,:)*0)*Ctr';%/h;
    DYLC=(y(ini+max(hz):end,:)-y(ini-h+max(hz):end-h,:)*0)*Ctr';%/h;
    DYfcastDYLC=DYfcastLC-DYLC;
        
    MSFELC(h,:)=mean((DYfcastDYLC).^2);
        

end
 
%Plotting the MSFE for the variables
figure(3);
cont = 1;
for ii=[1 2 3 4 6 7]
    subplot(3,3,cont); plot(MSFE(:,ii),'color',color,'LineStyle',LS,'LineWidth',LW); hold on;
    set(gca,'FontSize',12);
    title(variables{ii})
    if ismember(cont,[1 4 7]), ylabel('MSFE','FontSize',12); end
    if cont >= 7, xlabel('Quarters Ahead','FontSize',12); end
    cont = cont+1;
end
fig3=figure(3);
fig3.Position=[0 300   875   600];
% legend('PLR baseline','PLR invariant','PLR invariant (except C-Y)')
 
% Plotting the MSFE for the transformed variables
for ii=[2 3 5]
    subplot(3,3,cont); plot(MSFELC(:,ii),'color',color,'LineStyle',LS,'LineWidth',LW); hold on;
    set(gca,'FontSize',12);
    title(transformations{ii})
    if ismember(cont,[1 4 7]), ylabel('MSFE','FontSize',12); end
    if cont >= 7, xlabel('Quarters Ahead','FontSize',12); end
    cont = cont+1; 
end

% legend('PLR baseline','PLR invariant','PLR invariant (except C-Y)')
