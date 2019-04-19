%% This script will calculate 10 measures for all files containing RS in the folder: filepathName
% All files fcn*.m contain Matlab functions used in calculating measures.

%% This will suppress al Matlab warnings
warning('off','all')

%% Add path to use EEGLAB Matlab functions; Change path to your local copy of EEGLab
addpath(genpath('../eeglab14_1_2b/'));

%% Flag indicating number of channels for processing
% If flag1020 = 1 then we process only 10/20 channels according to p. 7 in HydroCelGSN_10-10.pdf
% If flag1020 = 0 then we process all channels accordingly.
flag1020 = 1;  

%% Flag indicating whether to use FIR
% If flagFiltered = 1 then FIR is applied.
% If flagFiltered = 0 then FIR is not used.
flagFiltered = 0; 

%% Downsample rate only samples every x values to reduce computation time for testing. Make '1' for max.
if exist('downsampleRate', 'var')
    disp(['downsampleRate = ', num2str(downsampleRate)])
else
    downsampleRate = 5; 
    disp(['downsampleRate not set... will use ', num2str(downsampleRate)])
end

%% Get file(s)
myFolderInfo = dir('../AllRAWfiles/Pilots/*3rs.RAW'); 
myFolderInfo = myFolderInfo(~cellfun('isempty', {myFolderInfo.date}));

% time stats
time_CD_PK = 0.;
time_MFDFA = 0.;
time_PSVG = 0.;
time_tot = tic;

%% Iterate through available files in the folder
for iFile = 1:size(myFolderInfo,1)
    disp([' File: ', num2str(iFile), ' ', myFolderInfo(iFile).name])   % File for processing
    
    %% Read binary simple Netstation file
    filename = myFolderInfo(iFile).name; 
    EEG = pop_readegi(['../AllRAWfiles/Pilots/', filename], [],[],'auto');

    %% Correct delay 
    EEG = correctDelay(EEG,22);

    %% Edit channel locations 
    myChanLocs = 'GSN-HydroCel-129.sfp';
    EEG = pop_chanedit(EEG, 'load',{myChanLocs 'filetype' 'autodetect'},'setref',{'4:128' 'Cz'},'changefield',{132 'datachan' 0});

    %% Re-reference and add 'Cz' back to the data
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'Cz'},'Y',{0},'X',{5.4492e-16},'Z',{8.8992},'sph_theta',{0},'sph_phi',{90},'sph_radius',{8.8992},'theta',{0},'radius',{0},'type',{''},'ref',{'Cz'},'urchan',{132},'datachan',{0}));

    %% Filter the data; 0.1 for low and 50 for high
    if flagFiltered==1
        EEG = pop_eegfiltnew(EEG, 0.1,50,33000,0,[],1);    
    end
    
    %% Correct DINs
    EEG.event = cleanTriggers_v3(EEG.event);

    %% Plot for checking
    %pop_eegplot( EEG, 1, 1, 1);

    %% Run ICA
    %EEG = pop_runica(EEG, 'extended',1,'interupt','on');

    %% Use for checking consistency of dataset
    EEG = eeg_checkset(EEG);

    %% Save dataset; 
    %EEG = pop_saveset( EEG, 'filename',strrep(filename,'.RAW','RS.set'),'filepath',filepathName);

    %% Calculate fractal dimensions and save the output to Excel spreadsheet
    % Prepare table for output; allocate memory
    tableOutput = struct2table(EEG.event);
    for jChan=1:size(EEG.chanlocs,2)
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_CD')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_PK')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_FNNB')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_D')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_LE')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_HFD')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_MSE')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_MFDFA_DQFIRST')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_MFDFA_MAXDQ')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_MFDFA_DQLAST')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_MFDFA_MAXMIN')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_LZ')) = {0}; 
        tableOutput(1, strcat(EEG.chanlocs(jChan).labels,'_VG')) = {0}; 
    end
    %% We process EOEC only 
    % The following epochs are used:
    % [DIN0 DIN0], [DIN1 - 30000 DIN1 + 30000], [DIN0 DIN1] and [DIN1 DIN0]
    
    % EOEC; 
    eventVec = 1:3;

    %% Iterate through events
    for iEvent = eventVec
         % Extract epochs
            % EOEC
            if size(tableOutput,1)<=5 
                switch iEvent
                    case 1 % [DIN1 - 60000 DIN1 + 60000],
                         tempDataAll = EEG.data(:, EEG.event(2).latency - 60000:EEG.event(2).latency + 60000);
                    case 2 % [DIN0 DIN1]
                         tempDataAll = EEG.data(:, EEG.event(1).latency:EEG.event(2).latency);
                    case 3 % [DIN1 DIN0]
                         tempDataAll = EEG.data(:, EEG.event(2).latency:EEG.event(3).latency);
                end
            end

            % Store results in this matrix for parallel processing purposes
            resultMat = zeros(size(EEG.chanlocs,2),13);
            
            % Select channels accroding to flag1020
            channelVec = []; % Initiate the variable
            if flag1020 == 1
                channelVec = [36, 104, 129, 24, 124, 33, 122, 22, 9, 14, 21, ...
                    15, 11, 70, 83, 52, 92, 58, 96, 45, 108];
            else
                channelVec = 1:size(EEG.chanlocs,2);
            end
            
        for jChan = 1:size(EEG.chanlocs,2);
            tic;
            
            if sum(channelVec==jChan)==1
                % CD, PK, FNNB, D
                tic;
                uf = 1; % Use fnn
                tt = 0; % Measure time - tic toc
                prt = 0; % Print results
                [CD, PK, FNNB, D] = fcnEMBED_v3(downsample(tempDataAll(jChan,:),downsampleRate),uf,tt,prt); 
                time_CD_PK = time_CD_PK + toc;	

                % Lyapunov Spectrum
                LE = fcnLE(downsample(tempDataAll(jChan,:)',downsampleRate),1);
 
                % Higuchi FD
                kmax = 5;
                HFD = fcnHFD(downsample(tempDataAll(jChan,:),downsampleRate), kmax);

                % MSE
                [MSE, ~, ~] = fcnSE(downsample(tempDataAll(jChan,:),downsampleRate));

                % MFDFA
                scmin = 16;
                scmax = 1024;
                scres = 19;
                exponents = linspace(log2(scmin),log2(scmax),scres);
                scale = round(2.^exponents);
                q = linspace(-5,5,101);
                m = 1;
                
                tic;  
                [Hq,tq,hq,Dq,Fq] = fcnMFDFA(downsample(tempDataAll(jChan,:),downsampleRate),scale,q,m,0);
                time_MFDFA = time_MFDFA + toc;
                
                % Remove NaN and Inf in Dq and hq
                hq = hq(~isnan(hq) | ~isinf(hq));
                Dq = Dq(~isnan(Dq) | ~isinf(Dq));
                % Create vector for storing outputs
                MFDFA = zeros(1,4);
                if (~isempty(hq) && ~isempty(Dq))
                    MFDFA(1) = Dq(1);
                    MFDFA(2) = max(Dq);
                    MFDFA(3) = Dq(end);
                    MFDFA(4) = max(hq)-min(hq);
                end
                
                % LZ
                LZ = fcnLZ(downsample(tempDataAll(jChan,:),downsampleRate)...
                    >= median(downsample(tempDataAll(jChan,:),downsampleRate)));
                
                % PSVG
                %tic;
                %VG = fcnPSVG_v3(downsample(tempDataAll(jChan,:),downsampleRate)');
                %time_PSVG = time_PSVG + toc;
                
                % Store results 
                if FNNB==0
                    FNNB = 1.0000e-16; % To avoid deleting the column later in the code
                end
                VG = 0;
                resultMat(jChan,:) = [CD, PK, FNNB, D, LE, HFD, MSE, MFDFA, LZ, VG];
            end
            
            % Show progress
            if rem(jChan, 10)==1
                disp([' Channel: ', num2str(jChan)])
            end
            toc
        end
        
        % Save output to a table
        resultMatT = resultMat';
        tableOutput{iEvent, 4:4 + 129*13 - 1} = resultMatT(:)'; 

        % Perform a checksum and display
        disp(['check nansum: ', num2str(nansum(resultMat(:)))])

	% Time stats
	disp([' CD_PK_v2 timing: ', num2str(time_CD_PK),...
	      ' MFDFA timing: ', num2str(time_MFDFA),...
	      ' total time: ', num2str(toc(time_tot)),...
	      ' [secs]'])
        
        % Save length of epoch
        tableOutput(iEvent, 'Epoch_length') = {size(tempDataAll,2)}; 
        disp([' Event: ', num2str(iEvent)])
    end % loop for events

    %% Save output to Excel spreadsheet for checking and futher processing
    % Leave only columns which have at least one non-zero value
    tableOutput = tableOutput(:,4:end);
    tableOutput = tableOutput(:,~all(tableOutput{:,:}==0));
    
    writetable(tableOutput,strrep(['../AllRAWfiles/Pilots/', filename],'.RAW','RS.xlsx'),'Sheet',1,'Range','A1')

end % loop for files



