% This script creates the .mat file for input based on a csv and some internal settings.  It is
% to allow replication without access to the internal databases.
%
% Kurt Lewis
% Created: December 3, 2018
%

% Clear the decks
clear all
clc

% Load the data in the csv file
[allData,hdr] = xlsread('../data/input/lvg.data.xls');

% We have loaded in the full data file, and we need to get the information that will
% allow us to initialize the state variables.  This information is essentially the
% trend component of the HP filter of the lgdp data starting from the beginning of the
% estimation sample

startInd = datefind(datenum(1961,03,31),x2mdate(allData(:,1)),1);

% The state initialization will begin using the calendar year before the estimation
% starts, so I need to get the trend component from the HP filter for 1960Q2-Q4, along
% with the growth rate of that trend (as measured by the first difference since this is
% log data).  Before I start, I need to get the trend and cycle from the full data
% beginning in 1960.
hpStart = startInd - 4;
gdpHP = allData(hpStart:end,2)*100;

% Get the trend and the cycle using quarterly settings for the HP filter, which is a
% smoothing parameter of 1600.  This follows the setup of HLW (2017). 
[Tr,Cy] = hpfilter(gdpHP,36000);
gHP = Tr(2:end) - Tr(1:end-1);
yPinit = Tr(4:-1:2);
gInit = gHP(3:-1:2);
zInit = zeros(2,1);

% Still use the same Initial conditions that HLW (2017) use, just for consistency
m0hp = [yPinit; gInit; zInit];

V0hp = [0.7 0.2 0.0 0.2 0.2 0.0 0.0;...
        0.2 0.2 0.0 0.0 0.0 0.0 0.0;...
        0.0 0.0 0.2 0.0 0.0 0.0 0.0;...
        0.2 0.0 0.0 0.2 0.2 0.0 0.0;...
        0.2 0.0 0.0 0.2 0.2 0.0 0.0;...
        0.0 0.0 0.0 0.0 0.0 0.2 0.2;...
        0.0 0.0 0.0 0.0 0.0 0.2 0.2];

%%% This is all the data
dates = x2mdate(allData(startInd:end,1));
gdp = allData(startInd:end,2)*100;
infl = allData(startInd:end,3); 
gdpl1 = allData(startInd:end,4)*100;
gdpl2 = allData(startInd:end,5)*100;
infll1 = allData(startInd:end,6);
infll2 = allData(startInd:end,7);
infll3 = allData(startInd:end,8);
infll4 = allData(startInd:end,9);
rl1 = allData(startInd:end,10);
rl2 = allData(startInd:end,11);
aX = zeros(size(gdp)); % Blanks, previously used to look at other data

% If we are using the EBP data, it comes in here

% Organize the data into matrices
y = [gdp';infl'];
u = [ones(size(gdpl1'));gdpl1';gdpl2';rl1';rl2';infll1'; (1/3)*(infll2+infll3+infll4)'; ...
         aX'];

%% Get some features of the data sorted out:
[m,T] = size(y)
[n,T1] = size(u)

% Organize the data into the structure
data.y = y;
data.u = u;
data.dates = dates;
% This vintage is based on the vintage of the data in the provided lvg.data.xls file
data.vintage = '2018-10-26';
data.m0 = m0hp;
data.V0 = V0hp;

% Save the file in the data/input directory
save('../data/input/lvg.data.new.mat','data')
