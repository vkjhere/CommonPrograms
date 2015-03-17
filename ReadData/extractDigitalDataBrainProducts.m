% This function is used to extract the digital data from raw brain products
% data files (.eeg, .vhdr, .vmrk)

% Each data file is characterized by four parameters - subjectName, expDate,
% protocolName and gridType.

% We assume that the raw data is initially stored in
% folderSourceString\data\rawData\{subjectName}{expDate}\

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to make sure that digital codes from different recording systems
% are similar, we convert the digital codes to the Blackrock format.

% The GRF Protocol automatically extends the length of two digital -
% trialStart and trialEnd, to 2 ms, making sure that they are captured
% properly as long as the sampling frequency exceeds 1 kHz. Further, the
% codes corresponding to reward are also recorded. The program therefore
% does the following:

% 1. Finds out which codes correspond to reward on/off, TrialStart and TrialEnd
% 2. Changes these codes to the format used in Blackrock
% 3. Ignores the other digital codes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [digitalTimeStamps,digitalEvents]=extractDigitalDataBrainProducts(subjectName,expDate,protocolName,folderSourceString,gridType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = [subjectName expDate protocolName '.vhdr'];
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
makeDirectory(folderName);
folderIn = fullfile(folderSourceString,'data','rawData',[subjectName expDate]);
folderExtract = fullfile(folderName,'extractedData');
makeDirectory(folderExtract);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use EEGLAB plugin "bva-io" to read the file
eegData = pop_loadbv(folderIn,fileName,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Digital Codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[digitalTimeStamps,digitalEvents] = getStimulusEvents(eegData.event);
disp(['Ditial events: Total: ' num2str(length(eegData.event)) ', Stimulus: ' num2str(length(digitalTimeStamps))]);

% We only consider codes that are separated by at least deltaLimit ms to make sure
% that none of the codes are during the transition period.

digitalTimeStamps = digitalTimeStamps/eegData.srate; % Time in seconds

deltaLimit = 3/1000; % ms 
dt = diff(digitalTimeStamps);
badDTPos = find(dt<=deltaLimit);

if ~isempty(badDTPos)
    disp([num2str(length(badDTPos)) ' of ' num2str(length(digitalTimeStamps)) ' (' num2str(100*length(badDTPos)/length(digitalTimeStamps),2) '%) are separated by less than ' num2str(deltaLimit) ' s and will be discarded']);
    digitalTimeStamps(badDTPos)=[];
    digitalEvents(badDTPos)=[];
end
end

function [eventTimes,eventVals] = getStimulusEvents(allEvents)

% All 16 digital pins of BrainProducts should be set to 'Stimulus'. If that
% is not the case, throw an error for now.

count=1;
for i=1:length(allEvents)
    if strcmp(allEvents(i).code, 'Response')
        error('Response pins should be set to Stimulus in the Brain Products configurations settings.');
    end   
    if strcmp(allEvents(i).code, 'Stimulus')
        eventTimes(count) = allEvents(i).latency;
        eventVals(count)   = str2double(allEvents(i).type(2:end)); %#ok<*AGROW>
        count=count+1;
    end
end

% MSB is set to negative. Change to positive
x = find(eventVals<0);
eventVals(x) = 2^16 + eventVals(x);
end