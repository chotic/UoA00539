function eventCorrect = cleanTriggers_v3(event)
%cleanTriggers takes EEG.event input and returns a cleaned version:
% - replace sequence [DIN8, DIN4, DIN2, DIN1] with DIN0 for RMET and BCST
% - replace sequence [DIN8, DIN4, DIN2] with DIN0 for Eyes Open/Close
% - leave only one DIN0 at the begining and at the end of the events

%% Loops through all the events contained in event; RMET and BCST
% If the sequence [DIN8, DIN4, DIN2, DIN1] is present replace 
% with [DIN0, DIN0, DIN0, DIN0]
    
for i = 1:(size(event,2)-3)
    if (strcmp(event(i).type,'DIN8') && strcmp(event(i+1).type,'DIN4') && ...
            strcmp(event(i+2).type,'DIN2') && strcmp(event(i+3).type,'DIN1') && ... 
            event(i).latency == event(i+1).latency && event(i+1).latency == event(i+2).latency && ...
            event(i+2).latency == event(i+3).latency) 
        % Replace with DIN0
        event(i).type = 'DIN0';
        event(i+1).type = 'DIN0';
        event(i+2).type = 'DIN0';
        event(i+3).type = 'DIN0';
    end
    
    if (strcmp(event(i).type,'DIN1') && strcmp(event(i+1).type,'DIN2') && ...
            strcmp(event(i+2).type,'DIN4') && strcmp(event(i+3).type,'DIN8') && ...
            event(i).latency == event(i+1).latency && event(i+1).latency == event(i+2).latency && ...
            event(i+2).latency == event(i+3).latency) 
        % Replace with DIN0
        event(i).type = 'DIN0';
        event(i+1).type = 'DIN0';
        event(i+2).type = 'DIN0';
        event(i+3).type = 'DIN0';
    end
    
    
end

for i = 1:(size(event,2)-2)
    if (strcmp(event(i).type,'DIN1') && strcmp(event(i+1).type,'DIN2') && strcmp(event(i+2).type,'DIN4') && ...
            event(i).latency == event(i+1).latency && event(i+1).latency == event(i+2).latency)
        % Replace with DIN0
        event(i).type = 'DIN0';
        event(i+1).type = 'DIN0';
        event(i+2).type = 'DIN0';
    end
    if (strcmp(event(i).type,'DIN4') && strcmp(event(i+1).type,'DIN2') && strcmp(event(i+2).type,'DIN1') && ...
            event(i).latency == event(i+1).latency && event(i+1).latency == event(i+2).latency)
        % Replace with DIN0
        event(i).type = 'DIN0';
        event(i+1).type = 'DIN0';
        event(i+2).type = 'DIN0';
    end
end
    
%% Loops through all the events contained in event; Eyes Open -> Eyes Closed
% If the sequence [DIN8, DIN4, DIN2] is present replace 
% with [DIN0, DIN0, DIN0] - ending sequence

% Apply only if there are less than 10 DINs to make sure that we work with
% Eyes Open -> Eyes Closed files
if size(event,2) <=10
    % Starting DIN
    if strcmp(event(1).type,'DIN2')
        % Replace with DIN0
        event(1).type = 'DIN0';
    end

    % Endind sequence
    for i = 1:(size(event,2)-2)
        if (strcmp(event(i).type,'DIN8') && strcmp(event(i+1).type,'DIN4') && ...
                strcmp(event(i+2).type,'DIN2')) || (strcmp(event(i).type,'DIN2') && ...
                strcmp(event(i+1).type,'DIN4') && strcmp(event(i+2).type,'DIN8')) 
            % Replace with DIN0
            event(i).type = 'DIN0';
            event(i+1).type = 'DIN0';
            event(i+2).type = 'DIN0';
        end
    end
end


%% Replace first block of DIN0 with only one DIN0 with the largest latency
% Replace last block of DIN0 with only one DIN0 with the smallest latency
indDIN0 = arrayfun(@(x) strcmp(x.type, 'DIN0'), event);

% Get starting and endind index
tempStart = [1, find(diff(indDIN0)==-1)];
tempStart = tempStart(~isnan(tempStart));
startInd = max(tempStart);
tempEnd = [size(event,2), find(diff(indDIN0)==1) + 1];
tempEnd = tempEnd(~isnan(tempEnd));
endInd = min(tempEnd);

% Create the correct version of event
ind = 1;
for i=1:endInd
    if (i >= startInd && i <= endInd)
        eventCorrect(ind) = event(i);
        ind = ind + 1;
    end
end   

end