% In this script we perform EDA and statistical analysis of RS, BCST and
% RMET files.


%% Analysis of files from Unprocessed_RS_BCST_RMET folder
computerName = {'pc', 'mac'};
measureName = {'CD', 'PK', 'LE', 'HFD', 'MSE_1', 'MSE_2', 'MSE_3', 'MSE_4',...
    'MSE_5', 'MSE_6', 'MSE_7', 'MSE_8', 'MSE_9', 'MSE_10', 'MSE_11', 'MSE_12',...
    'MSE_13', 'MSE_14', 'MSE_15', 'MSE_16', 'MSE_17', 'MSE_18', 'MSE_19', 'MSE_20',...
    'MFDFA_DQFIRST', 'MFDFA_MAXDQ', 'MFDFA_DQLAST', 'MFDFA_MAXMIN', 'KC'};

% Iterate through computerName and then measures
for iComp = 1:2
    % Iterate through files to collect data for plotting and analysis - RS, 
    % BCST, RMET
    fileList = dir(['../OutputFiles/Unprocessed_RS_BCST_RMET/',...
        computerName{iComp}, '*.xlsx']); 
    fileList = fileList(~cellfun('isempty', {fileList.date}));

    % Iterate through files to collect data for plotting and analysis
    for jFile = 1:length(fileList(:))
        filename = fileList(jFile).name; 
        filenameSplit = strsplit(filename, '_');
        fileTask = filenameSplit{3};
        
        % Read file and add computerName, filename and resultType to table
        fileTable = readtable(['../OutputFiles/Unprocessed_RS_BCST_RMET/',...
            filename]);
        % Delete rows with zeros only
        fileTable(~any(fileTable{:,:}, 2), :) = [];  %rows
             
        %  Choose type of processing based on file
        switch fileTask(end-1:end)
            case 'rs' % rs short 5 seconds
                fileTable([1, 3],:) = []; %Save only EO
                fileTable.Event = {'EO'};
                 % Add properties to data       
                fileTable.Type = repmat('Unprocessed', 1, 1);
                fileTable.Filename = {filenameSplit{3}};
                fileTable.Computer = repmat(computerName{iComp}, 1, 1);
            case 'et' % rmet
                % Take median across columns
                fileTableMedian = median(fileTable{:, :});
                
                % Replace data with mean and median
                fileTable(1:end,:) = [];
                fileTable{1,:} = fileTableMedian;
                
                % Add properties to data       
                fileTable.Type = repmat('Unprocessed', 1, 1);
                fileTable.Filename = {filenameSplit{3}};
                fileTable.Computer = repmat(computerName{iComp}, 1, 1);
                fileTable.Event = {'MEDIAN_RMET'};
            case 'st' % bcst  
                % Take median across columns
                fileTableMedian = median(fileTable{:, :});
                
                % Replace data with mean and median
                fileTable(1:end,:) = [];
                fileTable{1,:} = fileTableMedian;
                
                % Add properties to data       
                fileTable.Type = repmat('Unprocessed', 1, 1);
                fileTable.Filename = {filenameSplit{3}};
                fileTable.Computer = repmat(computerName{iComp}, 1, 1);
                fileTable.Event = {'MEDIAN_BCST'};
        end
        
        % Concatenate results
        if jFile == 1
            resultTable  = fileTable;
        else
            resultTable = [resultTable; fileTable];
        end
        
    end
    
    % Iterate through measures
    for jMeasure = 1:size(measureName, 2)
        % Plot boxplot for each measure
        indMeasure = find(~cellfun(@isempty, ...
            strfind(resultTable.Properties.VariableNames, measureName{jMeasure})));

        % Create example data
        EO = resultTable{find(strcmp(resultTable.Event, 'EO')), indMeasure}';
        RMET = resultTable{find(strcmp(resultTable.Event, 'MEDIAN_RMET')), indMeasure}';
        BCST = resultTable{find(strcmp(resultTable.Event, 'MEDIAN_BCST')), indMeasure}';
                
        % Prepare data for plotting
        data = cell(3, 3);
        for i = 1:3 %size(data, 1)
            PP1{i} = EO(:, i);
            PP2{i} = RMET(:, i);
            PP3{i} = BCST(:, i);
        end
               
        data = vertcat(PP1, PP2, PP3);

        % Perform tests
        for i = 1:3
            % Anderson-Darling goodness-of-fit hypothesis test 
            [hAD(i), pAD(i)] = adtest(EO(:, i) - RMET(:, i));

            % Returns the p-value of a paired, two-sided test for the null hypothesis 
            % that x – y comes from a distribution with zero median. 
            % Wilcoxon signed rank test
            [pWSR(i), hWSR(i), ~] = signrank(EO(:, i), RMET(:, i), 'tail', 'right', 'method', 'exact');

            % Returns a test decision for the paired-sample t-test
            [hTT(i), pTT(i), ~, ~] = ttest(EO(:, i), RMET(:, i), 'Tail', 'right');
        end

        % Plot boxplots with results
        pAD = round(pAD, 2); pWSR = round(pWSR, 2); pTT = round(pTT, 2); % round
        xlab = {['1003rs, ', num2str(pAD(1)), ', ', num2str(pWSR(1)), ', ', num2str(pTT(1))],...
            ['2003rs, ', num2str(pAD(2)), ', ', num2str(pWSR(2)), ', ', num2str(pTT(2))],...
            ['3003rs, ', num2str(pAD(3)), ', ', num2str(pWSR(3)), ', ', num2str(pTT(3))]};
        col = [102, 255, 255, 200; 
              51, 153, 255, 200;
              0, 0, 255, 200]/255;

        figure(1)
        fcnBoxplot(data', xlab, {'EO', 'MEDIAN\_RMET', 'MEDIAN\_BCST'}, col')
        title(strrep(['Unprocessed_RS_BCST_RMET, ', computerName{iComp}, ', ',...
            measureName{jMeasure}], '_', '\_'))

        % Save results
        saveas(gcf, ['../Results/Unprocessed_RS_BCST_RMET_', computerName{iComp}, '_', ...
            measureName{jMeasure}, '.png'])
    end %jMeasure 
end %iComp


