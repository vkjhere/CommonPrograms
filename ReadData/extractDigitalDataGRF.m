% This function is used to extract the digital data for the GRF Protocol.

% These are the following modes in which the GaborRFMap protocol has been used so far.

% 1. Target and mapping stimulus 0 are on and are presented synchronously,
% while mapping stimulus 1 is off. Digital codes are sent only for map0 stimulus.

% 2. The task is run in the fixation mode, in which target stimulus is off
% and all trials are catch trials. Digital codes are sent only for map0
% stimulus. The target is assumed to be synchronous with the mapping
% stimulus in this case.

function [goodStimNums,goodStimTimes,side] = extractDigitalDataGRF(folderExtract,ignoreTargetStimFlag,frameRate)

if ~exist('ignoreTargetStimFlag','var');   ignoreTargetStimFlag=0;      end
if ~exist('frameRate','var');              frameRate=100;               end

stimResults = readDigitalCodesGRF(folderExtract,frameRate); % writes stimResults and trialResults
side = stimResults.side;
[goodStimNums,goodStimTimes] = getGoodStimNumsGRF(folderExtract,ignoreTargetStimFlag); % Good stimuli
save(fullfile(folderExtract,'goodStimNums.mat'),'goodStimNums','goodStimTimes');
end

% GRF Specific protocols
function [stimResults,trialResults,trialEvents] = readDigitalCodesGRF(folderOut,frameRate)

if ~exist('frameRate','var');              frameRate=100;               end
kForceQuit=7;
digitalCodeLoss = 0; % assume no loss of digital codes

% Get the values of the following trial events for comparison with the dat
% file from lablib
trialEvents{1} = 'TS'; % Trial start
trialEvents{2} = 'TE'; % Trial End

load(fullfile(folderOut,'digitalEvents.mat'));

allDigitalCodesInDec = [digitalCodeInfo.codeNumber];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the times and values of the events in trialEvents

for i=1:length(trialEvents)
    pos = find(convertStrCodeToDec(trialEvents{i})==allDigitalCodesInDec);
    if isempty(pos)
        warning(['Code ' trialEvents{i} ' not found!!']);
    else
        trialResults(i).times = [digitalCodeInfo(pos).time]; %#ok<*AGROW>
        trialResults(i).value = [digitalCodeInfo(pos).value];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stimulus properties
azimuth          = [digitalCodeInfo(find(convertStrCodeToDec('AZ')==allDigitalCodesInDec)).value];
elevation        = [digitalCodeInfo(find(convertStrCodeToDec('EL')==allDigitalCodesInDec)).value];
contrast         = [digitalCodeInfo(find(convertStrCodeToDec('CO')==allDigitalCodesInDec)).value];
temporalFrequency= [digitalCodeInfo(find(convertStrCodeToDec('TF')==allDigitalCodesInDec)).value];
radius           = [digitalCodeInfo(find(convertStrCodeToDec('RA')==allDigitalCodesInDec)).value];
sigma            = [digitalCodeInfo(find(convertStrCodeToDec('SI')==allDigitalCodesInDec)).value];
spatialFrequency = [digitalCodeInfo(find(convertStrCodeToDec('SF')==allDigitalCodesInDec)).value];
orientation      = [digitalCodeInfo(find(convertStrCodeToDec('OR')==allDigitalCodesInDec)).value];

% Get timing
trialStartTimes    = [digitalCodeInfo(find(convertStrCodeToDec('TS')==allDigitalCodesInDec)).time];
taskGaborTimes     = [digitalCodeInfo(find(convertStrCodeToDec('TG')==allDigitalCodesInDec)).time];
mapping0Times      = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).time];
mapping1Times      = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).time];
mappingPlaidTimes  = [digitalCodeInfo(find(convertStrCodeToDec('MP')==allDigitalCodesInDec)).time];
numTrials = length(trialStartTimes);

% Get values
taskGaborValues    = [digitalCodeInfo(find(convertStrCodeToDec('TG')==allDigitalCodesInDec)).value];
mapping0Values     = [digitalCodeInfo(find(convertStrCodeToDec('M0')==allDigitalCodesInDec)).value];
mapping1Values     = [digitalCodeInfo(find(convertStrCodeToDec('M1')==allDigitalCodesInDec)).value];
mappingPlaidValues = [digitalCodeInfo(find(convertStrCodeToDec('MP')==allDigitalCodesInDec)).value];


if isempty(azimuth) || isempty(elevation) || isempty(contrast) || isempty(temporalFrequency) ...
        || isempty(radius) || isempty(sigma) || isempty(spatialFrequency) || isempty(orientation)
    
    disp('Digital codes for the stimuli are not sent. Read from Lablib data file later ...');
else
    
    % Check the default case - only mapping0/1 is on, and only its stimulus properties are put out.
    
    if (max(abs(diff([length(azimuth) length(elevation) length(contrast) length(temporalFrequency) ...
            length(radius) length(sigma) length(spatialFrequency) length(orientation)]))) > 0 )
        % This could happen when Blackrock crashed right in the middle of
        % the digital codes being sent for a stimulus
        warning('Length of stimulus properties are not even! Digital codes may have been lost, trying to recover...');
        digitalCodeLoss = 1;
        numStimProps = min([length(azimuth) length(elevation) length(contrast) length(temporalFrequency) ...
            length(radius) length(sigma) length(spatialFrequency) length(orientation)]);
        
        if abs(numStimProps - length(mapping0Times)) == 1 && isempty(mapping1Times)
            disp('Only Mapping 0 is used');
            stimResults.azimuth = convertUnits(azimuth(1:numStimProps)',100);
            stimResults.elevation = convertUnits(elevation(1:numStimProps)',100);
            stimResults.contrast = convertUnits(contrast(1:numStimProps)',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency(1:numStimProps)',10);
            stimResults.radius = convertUnits(radius(1:numStimProps)',100);
            stimResults.sigma = convertUnits(sigma(1:numStimProps)',100);
            stimResults.orientation = convertUnits(orientation(1:numStimProps)');
            stimResults.spatialFrequency = convertUnits(spatialFrequency(1:numStimProps)',100);

        elseif abs(numStimProps - length(mapping1Times)) == 1 && isempty(mapping0Times)
            disp('Only Mapping 1 is used');
            stimResults.azimuth = convertUnits(azimuth(1:numStimProps)',100);
            stimResults.elevation = convertUnits(elevation(1:numStimProps)',100);
            stimResults.contrast = convertUnits(contrast(1:numStimProps)',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency(1:numStimProps)',10);
            stimResults.radius = convertUnits(radius(1:numStimProps)',100);
            stimResults.sigma = convertUnits(sigma(1:numStimProps)',100);
            stimResults.orientation = convertUnits(orientation(1:numStimProps)');
            stimResults.spatialFrequency = convertUnits(spatialFrequency(1:numStimProps)',100);

        else
            disp('Both Mapping 0 and 1 are used');
            numStimProps = floor(numStimProps/2)*2; % should be even, since there are paired properties for Plaid
            stimResults.azimuth = convertUnits(azimuth(1:numStimProps)',100);
            stimResults.elevation = convertUnits(elevation(1:numStimProps)',100);
            stimResults.contrast = convertUnits(contrast(1:numStimProps)',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency(1:numStimProps)',10);
            stimResults.radius = convertUnits(radius(1:numStimProps)',100);
            stimResults.sigma = convertUnits(sigma(1:numStimProps)',100);
            stimResults.orientation = convertUnits(orientation(1:numStimProps)');
            stimResults.spatialFrequency = convertUnits(spatialFrequency(1:numStimProps)',100);
        end
        
    else

        if ((length(azimuth) == length(mapping0Times)) && isempty(mapping1Times))
            disp('Only Mapping 0 is used');
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);

        elseif ((length(azimuth) == length(mapping1Times)) && isempty(mapping0Times))
            disp('Only Mapping 1 is used');
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);

        else
            disp('Both Mapping 0 and 1 are used');
            stimResults.azimuth = convertUnits(azimuth',100);
            stimResults.elevation = convertUnits(elevation',100);
            stimResults.contrast = convertUnits(contrast',10);
            stimResults.temporalFrequency = convertUnits(temporalFrequency',10);
            stimResults.radius = convertUnits(radius',100);
            stimResults.sigma = convertUnits(sigma',100);
            stimResults.orientation = convertUnits(orientation');
            stimResults.spatialFrequency = convertUnits(spatialFrequency',100);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instruction trials
instructionTrials = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('IT')==allDigitalCodesInDec)).value])';
if length(instructionTrials) ~= numTrials 
    error('Number of instruction trial entries different from numTrials');
end

% Catch trials
catchTrials = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('CT')==allDigitalCodesInDec)).value])';
if length(catchTrials) ~= numTrials 
    error('Number of catch trial entries different from numTrials');
end

% TrialCertify & TrialEnd (eotCode)
% These two entries may be repeated twice during force quit
trialCertify = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TC')==allDigitalCodesInDec)).value])';
eotCodes = convertUnits([digitalCodeInfo(find(convertStrCodeToDec('TE')==allDigitalCodesInDec)).value])';

forceQuits = find(eotCodes==kForceQuit);
numForceQuits = length(forceQuits);

if length(eotCodes)-numForceQuits == numTrials
    disp(['numTrials: ' num2str(numTrials) ' numEotCodes: '  ...
        num2str(length(eotCodes)) ', ForceQuits: ' num2str(numForceQuits)]);
    goodEOTPos = find(eotCodes ~=kForceQuit);
    eotCodes = eotCodes(goodEOTPos);
    trialCertify = trialCertify(goodEOTPos);
else
    disp(['numTrials: ' num2str(numTrials) ' numEotCodes: '  ...
        num2str(length(eotCodes)) ', forcequits: ' num2str(numForceQuits)]);
    if numForceQuits > 0
        error('ForceQuit pressed after trial started'); % TODO - deal with this case
    else
        % For now, only extract trials which have an eotCode. Later compare
        % the trial start and end times with Lablib data.
        if numTrials > length(eotCodes)
            warning('More trialStart codes than trialEnd codes! Digital codes may have been lost, trying to recover...');
            % Limit number of trials till the last trial with an eotCode
            numTrials = length(eotCodes);
            lastValidTrialEndTime = digitalCodeInfo(find(convertStrCodeToDec('TE')==allDigitalCodesInDec)).time(end);
            numValidStims = length(find([digitalCodeInfo(find(convertStrCodeToDec('AZ')==allDigitalCodesInDec)).time] < lastValidTrialEndTime));
            % Check whether stimulus properties are also lost
            if digitalCodeLoss && numStimProps < numValidStims
                error(['Unrecoverable digital code loss! numStimProps: ' ...
                    num2str(numStimProps) ', numValidStims: ' num2str(numValidStims)]);
            end
            digitalCodeLoss = 1;
            % Trim trialResults
            trialStartTimes = trialStartTimes(1:numTrials);
            for i=1:length(trialEvents)
                trialResults(i).times = trialResults(i).times(1:numTrials);
                trialResults(i).value = trialResults(i).value(1:numTrials);
            end
            % Trim stimulus and timing
            stimResults.azimuth = stimResults.azimuth(1:numValidStims);
            stimResults.elevation = stimResults.elevation(1:numValidStims);
            stimResults.contrast = stimResults.contrast(1:numValidStims);
            stimResults.temporalFrequency = stimResults.temporalFrequency(1:numValidStims);
            stimResults.radius = stimResults.radius(1:numValidStims);
            stimResults.sigma = stimResults.sigma(1:numValidStims);
            stimResults.orientation = stimResults.orientation(1:numValidStims);
            stimResults.spatialFrequency = stimResults.spatialFrequency(1:numValidStims);            
            taskGaborTimes     = taskGaborTimes(taskGaborTimes < lastValidTrialEndTime);
            mapping0Times      = mapping0Times(mapping0Times < lastValidTrialEndTime);
            mapping1Times      = mapping1Times(mapping1Times < lastValidTrialEndTime);
            mappingPlaidTimes  = mappingPlaidTimes(mappingPlaidTimes < lastValidTrialEndTime);
            taskGaborValues    = taskGaborValues(taskGaborTimes < lastValidTrialEndTime);
            mapping0Values     = mapping0Values(mapping0Times < lastValidTrialEndTime);
            mapping1Values     = mapping1Values(mapping1Times < lastValidTrialEndTime);
            mappingPlaidValues = mappingPlaidValues(mappingPlaidTimes < lastValidTrialEndTime);
        else
            disp(['Something wrong with digital codes! numTrials: ' num2str(numTrials) ' numEotCodes: '  ...
                num2str(length(eotCodes)) ', ForceQuits: ' num2str(numForceQuits)]);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numStimTask  = getStimPosPerTrial(trialStartTimes, taskGaborTimes);
numStimMap0  = getStimPosPerTrial(trialStartTimes, mapping0Times);
numStimMap1  = getStimPosPerTrial(trialStartTimes, mapping1Times);
numStimMapPl = getStimPosPerTrial(trialStartTimes, mappingPlaidTimes);

% Check if Task and Mapping stim nums are the same for non-intruction trials
nonInstructionTrials = find(instructionTrials==0,numTrials);
% Handle Plaids case first: Map0 and Map1 stim nums should BOTH be
% the same as Task stim nums
if (max(abs(numStimTask(nonInstructionTrials) - numStimMap0(nonInstructionTrials)))==0) && ...
    (max(abs(numStimTask(nonInstructionTrials) - numStimMap1(nonInstructionTrials)))==0)
    disp('Mapping0/1 and Task times are the same');
    % Rudimentary checking to make sure we are dealing with Plaids
    if length(numStimMapPl(nonInstructionTrials)) ~= length(numStimMap0(nonInstructionTrials)) || ...
            length(numStimMapPl(nonInstructionTrials)) ~= length(numStimMap1(nonInstructionTrials))
        error('Number of Plaid stimuli does not match with Mapping0/1 stimuli');
    end
    if max(abs(numStimMapPl(nonInstructionTrials)/3 - numStimMap0(nonInstructionTrials))) == 0 && ...
            max(abs(numStimMapPl(nonInstructionTrials)/3 - numStimMap1(nonInstructionTrials))) == 0
        % Event mappingPlaid is generated once each for taskGabor, Map0 and
        % Map1 in the old code: the divisor 3 discounts for that
        numStims = numStimMapPl / 3;
        % Assume mappingPlaid stimuli time is the same as Map1 stimuli,
        % since Map1 comes on last to create the Plaid
        stimResults.time = mappingPlaidTimes(3:3:end)';
        % set flag to indicate old code
        oldCode = 1;
    elseif max(abs(numStimMapPl(nonInstructionTrials) - numStimMap0(nonInstructionTrials))) == 0 && ...
            max(abs(numStimMapPl(nonInstructionTrials) - numStimMap1(nonInstructionTrials))) == 0
        % Event mappingPlaid is generated once for every plaid stimulus
        numStims = numStimMapPl;
        stimResults.time = mappingPlaidTimes';
        oldCode = 0;
    else
        error('Something is wrong with Plaids protocol, please check!');
    end
    % Assign unique code for Plaid
    stimResults.side = 3;
    taskType = convertUnits(taskGaborValues)'; 
    
    if sum(taskType)==0 % Target is always null
        if oldCode
            taskType = convertUnits(mappingPlaidValues(3:3:end))';
        else
            taskType = convertUnits(mappingPlaidValues)';
        end
    end
    
elseif (max(abs(numStimTask(nonInstructionTrials) - numStimMap0(nonInstructionTrials)))==0)
    disp('Mapping0 and Task times are the same');
    numStims = numStimMap0;
    stimResults.time = mapping0Times';
    stimResults.side = 0;
    taskType = convertUnits(taskGaborValues)'; 
    
    if sum(taskType)==0 % Target is always null
        taskType = convertUnits(mapping0Values)';
    end
    
elseif (max(abs(numStimTask(nonInstructionTrials) - numStimMap1(nonInstructionTrials)))==0)
    disp('Mapping1 and Task times are the same');
    numStims = numStimMap1;
    stimResults.time = mapping1Times';
    stimResults.side = 1;
    taskType = convertUnits(taskGaborValues)';
    
    if sum(taskType)==0 % Target is always null
        taskType = convertUnits(mapping1Values)';
    end
    
else
    warning('Mapping0/1 and Task times not the same');
    
    if sum(numStimMap0)>0 && sum(numStimMap1)==0
        disp('Using Mapping0 times instead of task times...');
        numStimTask = numStimMap0;        % Assume task time is the same as mapping stimulus time                            
        numStims = numStimMap0;
        stimResults.time = mapping0Times';
        
        stimResults.side = 0;
        taskType = convertUnits(mapping0Values)'; % Assume task times and types are the same as M0
    
    elseif sum(numStimMap0)==0 && sum(numStimMap1)>0
        disp('Using Mapping1 times instead of task times...');
        numStimTask = numStimMap1;        % Assume task time is the same as mapping stimulus time
        numStims = numStimMap1;
        stimResults.time = mapping1Times';
        
        stimResults.side = 1;
        taskType = convertUnits(mapping1Values)'; % Assume task times and types are the same as M1
        
    else
        disp('Using MappingPlaid times instead of task times...');
        numStimTask = numStimMapPl;        % Assume task time is the same as mapping plaid time
        numStims = numStimMapPl;
        stimResults.time = mappingPlaidTimes';
        
        stimResults.side = 3;
        taskType = convertUnits(mappingPlaidValues)'; % Assume task times and types are the same as MP
        if max(abs(numStimMapPl(nonInstructionTrials)/3 - numStimMap0(nonInstructionTrials))) == 0 && ...
                max(abs(numStimMapPl(nonInstructionTrials)/3 - numStimMap1(nonInstructionTrials))) == 0
            % old code where event MP was generated thrice for each plaid stim
            numStimTask = numStimTask / 3;
            numStims = numStims / 3;
            stimResults.time = stimResults.time(3:3:end)';            
            taskType = taskType(3:3:end)';
        end
    end
end

% Okay... If there _was_ digital code loss, we have fudged around with the
% remaining codes. Do all remaining sanity checking here before saving
% information.
if digitalCodeLoss
    % check that we counted all stimuli properly
    if (stimResults.side < 2 && numValidStims ~= sum(numStims)) || ...
            (stimResults.side == 3 && numValidStims ~= 2*sum(numStims))
        error(['Unrecoverable digital code loss! numValidStims: ' ...
            num2str(numValidStims) ', sum(numStims): ' num2str(sum(numStims))]);
    end
end

posTask = 0;
pos=0;
for i=1:numTrials
    taskTypeThisTrial = taskType(posTask+1:posTask+numStimTask(i));
    if (numStims(i)>0)
        stimResults.type(pos+1:pos+numStims(i)) = taskTypeThisTrial;
        stimResults.trialNumber(pos+1:pos+numStims(i)) = i;
        stimResults.stimPosition(pos+1:pos+numStims(i)) = 1:numStims(i);

        if stimResults.side==0
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
                (mapping0Times(pos+1:pos+numStims(i)) - mapping0Times(pos+1))*frameRate;
        elseif stimResults.side==1
            stimResults.stimOnFrame(pos+1:pos+numStims(i)) = ...
                (mapping1Times(pos+1:pos+numStims(i)) - mapping1Times(pos+1))*frameRate;
        end

        stimResults.instructionTrials(pos+1:pos+numStims(i)) = instructionTrials(i); %always zero
        stimResults.catch(pos+1:pos+numStims(i)) = catchTrials(i);
        stimResults.eotCodes(pos+1:pos+numStims(i)) = eotCodes(i);
        stimResults.trialCertify(pos+1:pos+numStims(i)) = trialCertify(i);
        pos = pos+numStims(i);
    end
    posTask = posTask+numStimTask(i);
end

% Save in folderOut
save(fullfile(folderOut,'stimResults.mat'),'stimResults','digitalCodeLoss');
save(fullfile(folderOut,'trialResults.mat'),'trialEvents','trialResults');

end
function [goodStimNums,goodStimTimes] = getGoodStimNumsGRF(folderOut,ignoreTargetStimFlag)

if ~exist('ignoreTargetStimFlag','var');      ignoreTargetStimFlag=0;   end

load(fullfile(folderOut,'stimResults.mat'));

totalStims = length(stimResults.eotCodes);
disp(['Number of trials: ' num2str(max(stimResults.trialNumber))]);
disp(['Number of stimuli: ' num2str(totalStims)]);

% exclude uncertified trials, catch trials and instruction trials
tc = find(stimResults.trialCertify==1);
it = find(stimResults.instructionTrials==1);
ct = find(stimResults.catch==1);

if ignoreTargetStimFlag
    badStimNums = [it tc]; % catch trials are now considered good
else
    badStimNums = [it tc ct];
end

%eottypes
% 0 - correct, 1 - wrong, 2-failed, 3-broke, 4-ignored, 5-False
% Alarm/quit, 6 - distracted, 7 - force quit
%disp('Analysing correct, wrong and failed trials');
%badEOTs = find(stimResults.eotCodes>2); 
%disp('Analysing correct and wrong trials')
%badEOTs = find(stimResults.eotCodes>1); 
disp('Analysing only correct trials')
badEOTs = find(stimResults.eotCodes>0); 
badStimNums = [badStimNums badEOTs];

goodStimNums = setdiff(1:totalStims,unique(badStimNums));

% stim types
% 0 - Null, 1 - valid, 2 - target, 3 - frontpadding, 4 - backpadding
if ~ignoreTargetStimFlag
    disp('Only taking valid stims ');
    validStims = find(stimResults.type==1);
    goodStimNums = intersect(goodStimNums,validStims);
    
    %%%%%%%%%%%%%% Remove bad stimuli after target %%%%%%%%%%%%%%%%%%%%
    
    clear trialNums stimPos
    trialNums = stimResults.trialNumber(goodStimNums);
    stimPos   = stimResults.stimPosition(goodStimNums);
    
    % Get the target positions of the trialNums
    clear goodTrials
    goodTrials = unique(trialNums);
    
    clear targetPos
    for i=1:length(goodTrials)
        allStimWithThisTrialNum = find(stimResults.trialNumber==goodTrials(i));
        
        if sum(stimResults.catch(allStimWithThisTrialNum))>0        % catch trials
            targetPos(trialNums==goodTrials(i)) = inf; %#ok<*AGROW>
        else
            targetPos(trialNums==goodTrials(i)) = find(stimResults.type(allStimWithThisTrialNum)==2);
        end
    end
    
    validStimuliAfterTarget = find(stimPos>targetPos);
    if ~isempty(validStimuliAfterTarget)
        disp([num2str(length(validStimuliAfterTarget)) ' out of ' num2str(length(goodStimNums)) ' stimuli after target']);
        save(fullfile(folderOut,'validStimAfterTarget.mat'),'validStimuliAfterTarget');
    end
    
    goodStimNums(validStimuliAfterTarget)=[];
end
disp(['Number of good stimuli: ' num2str(length(goodStimNums))]);
goodStimTimes = stimResults.time(goodStimNums);
end
function outNum = convertUnits(num,f)

if ~exist('f','var');                       f=1;                        end

for i=1:length(num)
    if num(i) > 16384
        num(i)=num(i)-32768;
    end
end
outNum = num/f;
end
function [numStim,stimOnPos] = getStimPosPerTrial(trialStartTimes, stimStartTimes)

numTrials = length(trialStartTimes);

stimOnPos = cell(1,numTrials);
numStim   = zeros(1,numTrials);

for i=1:numTrials-1
    stimOnPos{i} = intersect(find(stimStartTimes>=trialStartTimes(i)),find(stimStartTimes<trialStartTimes(i+1)));
    numStim(i) = length(stimOnPos{i});
end
stimOnPos{numTrials} = find(stimStartTimes>=trialStartTimes(numTrials));
numStim(numTrials) = length(stimOnPos{numTrials});
end
