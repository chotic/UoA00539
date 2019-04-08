%% This script will evaluate measures for selected channel and time interval

%% This will suppress al Matlab warnings
warning('off','all')

%% Add path to use EEGLAB Matlab functions; Change path to your local copy of EEGLab
addpath(genpath('./'));

%% Compile mex code
mex computeDists.cpp
mex computeRatio.cpp
mex countGraphEdges.cpp

%% Get file you want to investigate
myFolderInfo = dir('../Raw files/Pilot3003.RAW'); 
myFolderInfo = myFolderInfo(~cellfun('isempty', {myFolderInfo.date}));
iFile = 1;

%% Read binary simple Netstation file
filename = myFolderInfo(iFile).name; 
EEG = pop_readegi(filename, [],[],'auto');

%% Correct delay 
EEG = correctDelay(EEG,22);

%% Edit channel locations 
myChanLocs = 'GSN-HydroCel-129.sfp';
EEG = pop_chanedit(EEG, 'load',{myChanLocs 'filetype' 'autodetect'},'setref',{'4:128' 'Cz'},'changefield',{132 'datachan' 0});

%% Re-reference and add 'Cz' back to the data
EEG = pop_reref( EEG, [],'refloc',struct('labels',{'Cz'},'Y',{0},'X',{5.4492e-16},'Z',{8.8992},'sph_theta',{0},'sph_phi',{90},'sph_radius',{8.8992},'theta',{0},'radius',{0},'type',{''},'ref',{'Cz'},'urchan',{132},'datachan',{0}));

%% Calculate measures
% Set global parameters
downsampleRate = 10; % Set downsampleRate
channel = 10; % Set channel for MFDFA calculations
epochStart = 10000; % Epoch start time
epochEnd = 100000; % Epoch end time

% Set parameters for MFDFA
scmin = 16;
scmax = 1024;
scres = 8;
exponents = linspace(log2(scmin),log2(scmax),scres);
scale = round(2.^exponents);
q = linspace(-5,5,101);
m = 1;
 
% Evaluate MFDFA for selected above parameters
[Hq,tq,hq,Dq,Fq] = fcnMFDFA(downsample(EEG.data(channel, epochStart:epochEnd),downsampleRate),scale,q,m,0);

% Extract Chris's parameters
MFDFA = zeros(1,4);
MFDFA(1) = Dq(1); % MFDFA_DQFIRST
MFDFA(2) = max(Dq); % MFDFA_MAXDQ
MFDFA(3) = Dq(end); % MFDFA_DQLAST
MFDFA(4) = max(hq) - min(hq); % MFDFA_MAXMIN

% Print Chris's parameters
MFDFA 

% Evaluate CD, PK and FNNB 
% Set parameters for CD, PK, FNNB
d = 10; % Maximum dimension
plt = 0; % Plot flag
tt = 1; % Tic flag
fnntol = 10; % Tolerance
slow = 0; % Slow flag 
usefnn = 1; % Use fnn flag
[CD, PK, FNNB] = fcnCD_PK_v2(downsample(EEG.data(channel, epochStart:epochEnd),downsampleRate),d,plt,tt,fnntol,slow,usefnn); 

% Print solution
CD
PK
FNNB

% Evaluate PSVG
VG = fcnPSVG(downsample(EEG.data(channel, epochStart:epochEnd),downsampleRate)');

% Print output
VG

%% Plot selected time interval
plot(EEG.data(channel, epochStart:epochEnd))