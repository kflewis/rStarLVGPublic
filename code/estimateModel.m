% To run a replication of the output in the files in data\output use this cede.  Otherwise leave
% it commented out.
jid = '10147376';
Wid = '1';
cede = str2double(strcat(jid, Wid));
rng(cede,'twister');

%% This is where outputs should be delivered.  Data outputs, things like the rStar and bounds estimates
outDataDir = '../data/output';

%% Switch for baseline vs alternative model in LVG (2018)
baseAlt = 1 % 0 = base, 1 = alt

%% MAT file output.
%% This cannot be a relative reference, it must be direct because the files generated
%% will be large and we can't risk dropping it into a home or production folder
%% on the server it is run on.
if ispc
    matOutDir = 'C:/scratch/';
else
    %matOutDir = '/scratch/';
    matOutDir = '/fst/scratch-m1kfl00/rStarPublic/';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This block controls whether it is the baseline or the alternative specification being
%% estimated 

%% File names and reporting information, Alternative Specification
if baseAlt == 1
    nameFile = ['AltModel_' datestr(now,'yyyymmdd') '_' datestr(now,'HHMMSS')];
    exten = '.mat';
    shortDes = 'AltModel'
    estLoc = [1 2 3 4 5 6 7 8 9 10 13]
elseif baseAlt == 0
    %%% File names and reporting information, Baseline Specification
    nameFile = ['BaseModel_' datestr(now,'yyyymmdd') '_' datestr(now,'HHMMSS')];
    exten = '.mat';
    shortDes = 'BaseModel'
    estLoc = [1 2 3 4 5 6 7 8 9 10]
else 
    error('Set which model is being estimated via baseAlt switch');
end

% The mat file that saves all the output
sf = [matOutDir,nameFile]

%% Run on a single machine, full size run
runInfo = ['500K burn in, 10 skip with 10,000 total draws']
burnIn = 500000
M = 10000
skip = 10

%%% Uncomment for shorter run, directory checking, diagnostics, etc.
%runInfo = ['250 burn in, 1 skip with 500 total draws, meant for TESTING']
%burnIn = 250
%M = 250
%skip = 1

totDraw = burnIn + (M*skip)

% The step size parameter for the MH algorithm when everything takes the same size step
% on average.  We will have random step sizes that average this size. 
stepSize = 0.02^2;

% How often (meaning after how many steps) to report information from the algorithm in the output
modShow = 1000;

% Bring in the vintaged data file
disp('Loading Data...');
load('../data/input/lvg.data.mat');
d1 = data;

% Identify the data structure of the model
y = d1.y;
u = d1.u;
dates = datenum(d1.dates);
m0 = d1.m0;
V0 = d1.V0;
vintage = d1.vintage;

%% Get some features of the data:
[m,T] = size(y);
[n,T1] = size(u);

%% Throw errors if the data isn't matching int erms of number of time periods
if (T~=T1)
    error('The number of observations for y and u are not the same')
end

% These are hypothetical values that can be set independently in case we don't want to estimate
% them.  The parameters being estimated are chosen using estLoc above.  These are approximately
% close to values found in the WP version of HLW.
a1pa2 = 0.95; % Alpha_1 + alpha_2
a2 = 0;       % alpha_2 unknown, set to zero for now 
a1 = a1pa2 - a2;
ar = -0.1;    
bY = 0.1;   
b1 = 0.7;     % b_pi 
s1 = 0.4;     % sig_ytilde in HLW
s2 = 0.8;     % sig_pi in HLW
s3 = 0.2;     % sig_z in HLW 
s4 = 0.6;     % sig_ystar in HLW
s5 = 0.12/4;  % sig_g in HLW
aX = 0;       % Extra exogenous variable potential coefficient can be set as needed

rhoG = 1;     % Picking 1 for both rhos, so that if we choose not to estimate them (removing estLoc 12 and 13)
rhoZ = 1;     % they are just the usual random walks.  The initial draw is determined by prior draw each time.

muZ = 0;      % Again, using 0 as the default value so that just not estimating it drops from the model
muG = 0;   % This is roughly the modal estimated value when both g and z are allowed to be
              % stationary is 0.65

% Approx. Lambda values from HLW
lambG = 0.05; 
lambZ = 0.03;

theta0 = [a1 a2 ar b1 bY s1 s2 s3 s4 s5 aX rhoG rhoZ muG muZ];

%% Keep the names of the variables
varNames = {'a1', 'a2', 'ar', 'b1', 'bY', 's1', 's2', 's3', 's4', 's5', 'aX', 'rhoG', ...
            'rhoZ', 'muG', 'muZ'};
varNamesLatex = {'$a_1$', '$a_2$', '$a_r$', '$b_1$', '$b_Y$', '$s_1$', '$s_2$',...
                 '$s_3$', '$s_4$', '$s_5$', '$a_X$', '$\rho_g$', ...
                 '$\rho_z$', '$\mu_g$', '$\mu_z$'};

% Construct the parameters based on the real info
[A0,B0,C0,D0,F0,G0,S0] = paramMat_distribute(theta0,y,u);

% Setting the parameters which will be estimated (estLoc) and those that will take assumed values (realLoc)    
realLoc = setdiff(1:15,estLoc);

% Draw from the prior
pDraw = priorDraw_distribute;

%% Initiate the collection of theta draws, and the step shock proerties
tDraws = NaN(totDraw,length(pDraw));
shkMu = zeros(size(pDraw));
shkSig = diag([1 1 1 1 1 .5 .5 500 .5 .25 1 1 100 1 1]);

% Except for the places we are NOT trying to estimate, here we put in the prior
theta0(estLoc) = pDraw(estLoc);

% Now, evaluate the prior for this theta
Lpi0 = -Inf;
while Lpi0 == -Inf
    Lpi0 = priorEval_distribute(theta0,estLoc);
    if Lpi0 == -Inf
        pDraw = priorDraw_distribute;
        theta0(estLoc) = pDraw(estLoc);
    end
end

% Construct the state space for the given parameters
[A,B,C,D,F,G,S] = paramMat_distribute(theta0,y,u);

% Based on the model, we need to get the size of the state variable
q = size(A,1);

%% Initiate the collection of the state variables and the likelihood value that corresponds
sDraws = NaN(q,T,totDraw);
mDraws = NaN(q,T,totDraw);
LLdraws = NaN(totDraw,1);
LpiDraws = NaN(totDraw,1);

% Evaluate the model with ffbs
[LL0,LLt,StDraw0,m,V,FFDraw0] = ffbs_kflfvg(y,u,A,B,C,D,F,G,S,m0,V0,1);

% Get the baseline R value
R0 = Lpi0 + LL0;

%% The first iteration of the MH loop
Rold = R0;
tOld = theta0;
sOld = StDraw0;
mOld = m;
VOld = V;
LLold = LL0;
LpiOld = Lpi0;
AR = NaN(totDraw,1);

LpiProp = -Inf;
while (LpiProp == -Inf)
    
    % Create the proposal, we aren't going far from the prior here because we know it
    % will be within the range where the priorEval won't return -Inf.
    tProp = tOld + mvnrnd(shkMu,1e-5*stepSize*shkSig);
    
    % Only estimate the parameters we ask to estimate (for debugging)
    tProp(realLoc) = theta0(realLoc);
    
    % Get the log-prior evaluation, it may be neg. infinity
    LpiProp = priorEval_distribute(tProp,estLoc);
    
end

% Construct the state space structures for the given parameters
[A,B,C,D,F,G,S] = paramMat_distribute(tProp,y,u);

% Get the loglikelihood of the proposed parameter vector
[LLProp,LLtProp,StDrawProp,mProp,VProp] = ffbs_kflfvg(y,u,A,B,C,D,F,G,S,m0,V0,1);

% Assemble the proposed statistic
Rprop = LpiProp + LLProp;
R = Rprop - Rold;

% Get the comparison
w = log(rand);

% Checking the comparison, put either the proposed value or the old value into the
% vector for the first entry
if (R >= w)
    % Assign the proposed values to the array
    tDraws(1,:) = tProp;
    
    % Set the old to the proprosed
    Rold = Rprop;
    tOld = tProp;
    sOld = StDrawProp;
    mOld = mProp;
    VOld = VProp;
    LLold = LLProp;
    LpiOld = LpiProp;
    
    AR(1,1) = 1;
    
else

    % Assign the old values to the array
    tDraws(1,:) = tOld;                     
    AR(1,1) = 0;

end

% Clear the temp variables
clear A B C D F G S 

% Now that we have the initial setup for the MH loop, we can run through the rest of
% the draws.
tStart = now();
tic
disp(sprintf('Entering Main Loop: %s',datestr(tStart)))

% Initialize skipping counter
j = 1;

for i = 2:totDraw
    
    LpiProp = -Inf;
    while (LpiProp == -Inf)
        
        % Create the proposal
        tProp = tOld + mvnrnd(shkMu,stepSize*shkSig);
        
        % Correct for the real stuff vs the estimated stuff 
        tProp(realLoc) = theta0(realLoc);                
        
        % Get the log-prior evaluation, it may be neg. infinity
        LpiProp = priorEval_distribute(tProp,estLoc);

    end

    % Construct the state space structures for the given parameters
    [A,B,C,D,F,G,S] = paramMat_distribute(tProp,y,u);
    
    % Get the loglikelihood of the proposed parameter vector
    [LLProp,LLtProp,StDrawProp,mProp,VProp] = ffbs_kflfvg(y,u,A,B,C,D,F,G,S,m0,V0,1);

    % Assemble the proposed statistic
    Rprop = LpiProp + LLProp;
    R = Rprop - Rold;

    % Get the comparison
    w = log(rand);
    
    % Checking the comparison, put either the proposed value or the old value into the
    % vector for the first entry
    if (R >= w)
                
        % Assign the proposed values to the array
        tDraws(i,:) = tProp;
   
        % Set the old to the proprosed
        Rold = Rprop;
        tOld = tProp;
        sOld = StDrawProp;
        mOld = mProp;
        VOld = VProp;
        LLold = LLProp;
        LpiOld = LpiProp;
        
        % Change the Accept-Reject
        AR(i,1) = 1;
    
    else

        % Assign the old values to the array
        tDraws(i,:) = tOld;                     
       
        % Change the Accept-Reject
        AR(i,1) = 0;
    end
    
    if (mod(i,modShow) == 0)

        rptStr = sprintf('Draw %i of %i, %3.1f percent.  Overall A/R: %3.1f, Last %i A/R: %3.1f',...
                          i,totDraw,(i/totDraw*100),(nanmean(AR)*100),modShow,...
                          (nanmean(AR((i-(modShow-1)):i,1))*100));
        
        disp(rptStr);
    end
    
    
    % Clear the temp variables
    clear A B C D F G S  
    clear LLProp LLtProp StDrawProp mProp VProp w R
    
end
tEnd = now();
disp(sprintf('Finished Main Loop: %s',datestr(tEnd)))
toc

%% At this point we want to extract the every skipth draw
tDrawsSkip = tDraws(burnIn+1:skip:end,:);

% Populate the res structure  ffbs_kflfvg(y,u,A,B,C,D,F,G,S,m0hp,V0hp,1);
res.dates = dates;
res.y = y;
res.u = u;
res.vintage = vintage;
res.m0 = m0;
res.V0 = V0;
res.tDraws = tDrawsSkip;
res.varNames = varNames;
res.varNamesLatex = varNamesLatex;

% Save all of the results in a mat file
sFileName = strcat(nameFile,exten);
sf = fullfile(matOutDir,sFileName);
save(sf,'-v7.3')

disp('Smaller results .mat file saved with RES structure.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 
%%%%% FROM HERE DOWN CAN BE RUN SEPERATELY STARTING FROM THE SMALLER RESULTS .MAT FILE TO
%%%%% GENERATE THE ESTIMATES OF THE UNOBSERVED VARIABLES.  THE SMALLER RESULTS .MAT FILE HAS THE
%%%%% DATA AND THE DRAWS FROM THE POSTERIOR, SO THEY CAN BE RUN THROUGH THE FILTER PROGRAM USING
%%%%% THIS CODE AND GENERATE THE UNOBSERVED COMPONENTS WITHOUT RUNNING THE ESTIMATION AGAIN.
%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we will always assume that we want both the filtered and the
% smoothed output
bs = 1;

% These are the vectors that extract the relevant state variables
% To extract potential output (pp), z (qqZ) and ANNUAL g (qqG) and full r*
pp = [1 0 0 0 0 0 0]'; 
qqZ = [0 0 0 0 0 1 0]'; 
qqG = [0 0 0 4 0 0 0]';
qq = [0 0 0 4 0 1 0]';

% The definition of long-term
lt = 40;

% Grab the results structure
r = res;

% Determine the important sizes
T = size(r.y,2);
Draws = size(r.tDraws,1);

% Add the data to the output
outVals.y = r.y;
outVals.u = r.u;
outVals.dates = r.dates;
outVals.vintage = r.vintage;
outVals.varNames = r.varNames;
outVals.varNamesLatex = r.varNamesLatex;
outVals.m0 = r.m0;
outVals.V0 = r.V0;

% Loop over the draws and build the estimates of the states
for i = 1:Draws
    
    if (mod(i,200)==0)
        disp(sprintf('Draw %i of %i, %2.1f Percent',i,Draws,(i/Draws*100)));
    end
    
    theta = r.tDraws(i,:);
    [A,B,C,D,F,G,S] = paramMat_distribute(theta,r.y,r.u);
    [LL,LLt,BSDraw,m,V,FFDraw] = ffbs_kflfvg(r.y,r.u,A,B,C,D,F,G,S,r.m0,r.V0,bs);
    outVals.theta(i,:) = theta;
    outVals.LL(i) = LL;
    outVals.LLt(:,i) = LLt;
    outVals.BSDraw(:,:,i) = BSDraw;
    outVals.FFDraw(:,:,i) = FFDraw;

    % The draws
    sTempBS = BSDraw;
    sTempFF = FFDraw;
    
    % Generate the gap from the mean state estimate
    outVals.yPBSDraws(:,i) = (pp'*sTempBS)';    
    outVals.yPFFDraws(:,i) = (pp'*sTempFF)';    
    outVals.zBSDraws(:,i) = (qqZ'*sTempBS)';    
    outVals.zFFDraws(:,i) = (qqZ'*sTempFF)';    
    outVals.gBSDraws(:,i) = (qqG'*sTempBS)';
    outVals.gFFDraws(:,i) = (qqG'*sTempFF)';
    outVals.rSBSDraws(:,i) = (qq'*sTempBS)';
    outVals.rSFFDraws(:,i) = (qq'*sTempFF)';
    
    % Need an additional loop to build the LR r* 
    muGdraws = theta(14);
    muZdraws = theta(15);
    rhoGdraws = theta(12);
    rhoZdraws = theta(13);
    
    for ii = 1:T
        outVals.rSLRBSDraws(ii,i) = muGdraws*(1-(rhoGdraws^(lt+1))) + ...
            (rhoGdraws^(lt+1))*outVals.gBSDraws(ii,i) + muZdraws*(1-(rhoZdraws^(lt+1))) + ...
            (rhoZdraws^(lt+1))*outVals.zBSDraws(ii,i);
        
        outVals.rSLRFFDraws(ii,i) = muGdraws*(1-(rhoGdraws^(lt+1))) + ...
            (rhoGdraws^(lt+1))*outVals.gFFDraws(ii,i) + muZdraws*(1-(rhoZdraws^(lt+1))) + ...
            (rhoZdraws^(lt+1))*outVals.zFFDraws(ii,i);
    end
    
    clear theta A B C D F G S LL LLt BSDraw m V FFDraw sTempBS sTempFF ...
        muGdraws muZdraws rhoGdraws rhoZdraws
end

%% Find the modes
for k = 1:size(outVals.theta,2)
    loc = k;
    [fs,xs] = ksdensity(outVals.theta(:,loc));
    [big,peakInd] = max(fs);
    pkSmooth(loc) = xs(peakInd);
    clear loc
end

%% All of the required draws for the lambdas
aRR = outVals.theta(:,3);
sOne = outVals.theta(:,6);
sTwo = outVals.theta(:,7);
sThree = outVals.theta(:,8);
sFour = outVals.theta(:,9);
sFive = outVals.theta(:,10);

lambGdraws = sFive./sFour;
lambZdraws = -aRR.*sThree./sOne;
a1pa2 = outVals.theta(:,1) + outVals.theta(:,2);

[fs,xs] = ksdensity(lambGdraws);
[big,peakInd] = max(fs);
pkSmooth(length(varNames)+1) = xs(peakInd);

[fs,xs] = ksdensity(lambZdraws);
[big,peakInd] = max(fs);
pkSmooth(length(varNames)+2) = xs(peakInd);

[fs,xs] = ksdensity(a1pa2);
[big,peakInd] = max(fs);
pkSmooth(length(varNames)+3) = xs(peakInd);

[fs,xs] = ksdensity(outVals.LL);
[big,peakInd] = max(fs);
pkSmooth(length(varNames)+4) = xs(peakInd);

% Add the lambdas to the table data
tab1Data = outVals.theta;
tabNames = varNames;

tab1Data = [tab1Data lambGdraws lambZdraws a1pa2 outVals.LL'];
tabNames(length(varNames)+1) = {'lambG'};
tabNames(length(varNames)+2) = {'lambZ'};
tabNames(length(varNames)+3) = {'a1+a2'};
tabNames(length(varNames)+4) = {'LogLike'};

estLoc2 = [estLoc (length(varNames)+1) (length(varNames)+2) (length(varNames)+3) (length(varNames)+4)];

% Spit out some properties of the parameter posterior distributions
for k = 1:size(tab1Data,2)
    loc = k;
    varX = tab1Data(:,loc);
    meanTab(loc) = mean(varX);
    medTab(loc) = median(varX);
    stdTab(loc) = std(varX);
    p5(loc) = prctile(varX,5);
    p10(loc) = prctile(varX,10);
    p50(loc) = prctile(varX,50);
    p90(loc) = prctile(varX,90);
    p95(loc) = prctile(varX,95);
    nameTab{loc} = tabNames(loc);
    
    clear loc varX
end

% To populate the LL values for the relevant parts of the table, we need the three values shown
% in the table
meanTheta = meanTab(1:size(outVals.theta,2));
medTheta = medTab(1:size(outVals.theta,2));
modeTheta = pkSmooth(1:size(outVals.theta,2));

[A,B,C,D,F,G,S] = paramMat_distribute(meanTheta,r.y,r.u);
[LLmean,LLt,StDraw,m,V] = ffbs_kflfvg(r.y,r.u,A,B,C,D,F,G,S,r.m0,r.V0,1);

[A,B,C,D,F,G,S] = paramMat_distribute(medTheta,r.y,r.u);
[LLmed,LLt,StDraw,m,V] = ffbs_kflfvg(r.y,r.u,A,B,C,D,F,G,S,r.m0,r.V0,1);

[A,B,C,D,F,G,S] = paramMat_distribute(modeTheta,r.y,r.u);
[LLmode,LLt,StDraw,m,V] = ffbs_kflfvg(r.y,r.u,A,B,C,D,F,G,S,r.m0,r.V0,1);

% Gather the LL values for the vectors given in the table
T1 = table(round(meanTab(estLoc2)',3),round(medTab(estLoc2)',3),...
         round(pkSmooth(estLoc2)',3),round(stdTab(estLoc2)',3),...
         round(p5(estLoc2)',3),round(p10(estLoc2)',3),round(p50(estLoc2)',3),...
         round(p90(estLoc2)',3),round(p95(estLoc2)',3),...
         'RowNames',tabNames(estLoc2)','VariableNames',...
         {'Mean','Median','Peak','StDev','five','ten','fifty','ninety','ninetyFive'});

llCell = num2cell([round([LLmean LLmed LLmode],2) NaN(1,6)]);
TLL = table(llCell{:},'RowNames',{'LogLike'},'VariableNames',{'Mean','Median','Peak','StDev','five','ten','fifty','ninety','ninetyFive'});
T1(end,:) = TLL

% Pull out the unobserved paths at percentiles ppp
ppp = [2.5 5 16 50 84 95 97.5];

outVals.rSBSpath = prctile(outVals.rSBSDraws,ppp,2);
outVals.rSFFpath = prctile(outVals.rSFFDraws,ppp,2);
outVals.rSLRBSpath = prctile(outVals.rSLRBSDraws,ppp,2);
outVals.rSLRFFpath = prctile(outVals.rSLRFFDraws,ppp,2);
outVals.yPBSpath = prctile(outVals.yPBSDraws,ppp,2);    
outVals.yPFFpath = prctile(outVals.yPFFDraws,ppp,2);    
outVals.zBSpath = prctile(outVals.zBSDraws,ppp,2);    
outVals.zFFpath = prctile(outVals.zFFDraws,ppp,2);    
outVals.gBSpath = prctile(outVals.gBSDraws,ppp,2);
outVals.gFFpath = prctile(outVals.gFFDraws,ppp,2);

%% Doing the log marginal data density calculation
LLmed = median(outVals.LL);
LMDD = LLmed - log(mean(exp(-outVals.LL+LLmed)));
disp(sprintf('Log Marginal Data Density Value: %3.1f \n', LMDD));

dateNames = strcat(num2str(year(outVals.dates)),'Q',num2str(quarter(outVals.dates)));
sh = 15;

rStarShow = outVals.rSBSpath((end-sh):end,4);
gShow = outVals.gBSpath((end-sh):end,4);
zShow = outVals.zBSpath((end-sh):end,4);
rStarLRShow = outVals.rSLRBSpath((end-sh):end,4);
dateShow = cellstr(dateNames((end-sh):end,:));

T2 = table(round(rStarShow,3),round(gShow,3),round(zShow,3),round(rStarLRShow,3),...
           'RowNames',dateShow,'VariableNames',{'rStar','g','z','LR_rStar'})

%% Generate the data files for smoothed and filtered output
dateFull = cellstr(dateNames);
baseN = {'rStar','rStarLR','yP','g','z'};
ind = 1;
for i = 1:length(baseN)
    for j = 1:length(ppp)
        colNames{ind} = strrep(strcat(baseN{i},'_',num2str(ppp(j))),'.','_');
        ind = ind + 1;
    end
end

tableDataFF = [outVals.rSFFpath outVals.rSLRFFpath outVals.yPFFpath outVals.gFFpath outVals.zFFpath];
tableDataBS = [outVals.rSBSpath outVals.rSLRBSpath outVals.yPBSpath outVals.gBSpath outVals.zBSpath];
        
TFF = array2table(round(tableDataFF,3),'RowNames',dateFull,'VariableNames',colNames);
TBS = array2table(round(tableDataBS,3),'RowNames',dateFull,'VariableNames',colNames);

writetable(TFF,'../data/output/lvgEstimatesFiltered.csv','WriteRowNames',true);
writetable(TBS,'../data/output/lvgEstimatesSmoothed.csv','WriteRowNames',true) 

%% Build the SDDR value
% Prior Centered at 1, sd = 0.5 truncated to be positive
rhoZprior = truncate(makedist('Normal','mu',1','sigma',0.5),0,Inf);
 
% Find the value of the pdf of the prior at 1
priorPdf = pdf(rhoZprior,1);

%% The marginal posterior density from rhoZ draws, using ksdensity.

% First, we are going to establish points for the ksdensity
ksDom = 0:0.001:1.2;

% Then we will build the posterior using ksdensity over those points
[fs,xs] = ksdensity(outVals.theta(:,13),ksDom);

% Then we find the posterior value at 1
postPdf = fs(1001);

% Build the BF
outVals.BF = priorPdf/postPdf;

sFileName = strcat(nameFile,'_outVals',exten);
sf2 = fullfile(matOutDir,sFileName);
save(sf2,'-v7.3')
disp('Larger results .mat file saved with OUTVALS structure.');

disp('Script complete.')
