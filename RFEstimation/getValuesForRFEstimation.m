% measure: 'LFP' or 'Spikes'
% timeRanges: a cell array of times - in seconds

function getValuesForRFEstimation(subjectName,expDate,protocolName,folderSourceString,gridType,measure,timeRanges,removeAvgRef,goodElectrodes)

if ~exist('timeRanges','var');          timeRanges=[];                   end
if ~exist('removeAvgRef','var');         removeAvgRef=0;                 end

stimPosGreaterThanOne=1;

if isempty(timeRanges)
    timeRanges{1}=[40 100]/1000;
    timeRanges{2}=[0 200]/1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% foldernames
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

folderOut1 = fullfile(folderName,'RFMeasures');
makeDirectory(folderOut1);
folderOut = fullfile(folderOut1,measure);
makeDirectory(folderOut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if removeAvgRef
    disp('Removing average reference');
    load(fullfile(folderSegment,'LFP','avgRef'));
    avgRef = analogData;
    fileTag = 'AvgRefRemoved';
else
    fileTag = '';
end

% Load stimulus parameters, bad trials
load(fullfile(folderExtract,'parameterCombinations.mat'));
aLength = length(aValsUnique);
eLength = length(eValsUnique);

numTimePeriods = length(timeRanges);

load(fullfile(folderSegment,'badTrials.mat'));

if strcmp(measure,'LFP')    % Case 1 - LFP analysis
    % Get Time Ranges
    load(fullfile(folderSegment,'LFP','lfpInfo.mat'));
    timePos = cell(1,numTimePeriods);
    for i=1:numTimePeriods
        timePos{i} = intersect(find(timeVals>=timeRanges{i}(1)),find(timeVals<timeRanges{i}(2)));
    end
    
    stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosGreaterThanOne);
    
    for i=1:length(analogChannelsStored)
        channelNumber = analogChannelsStored(i);
        
        % Get LFP data
        clear signal analogData meanLFPData numStimuli
        load(fullfile(folderSegment,'LFP',['elec' num2str(channelNumber) '.mat']));
        if removeAvgRef
            analogData = analogData-avgRef;
        end
        
        for a=1:aLength
            for e=1:eLength
                
                clear goodPos
                goodPos = intersect(parameterCombinations{a,e,1,end,end},stimPos); %#ok<*USENS>
                goodPos = setdiff(goodPos,badTrials);
                
                if isempty(goodPos)
                    rfValsRMS(e,a,channelNumber,1:numTimePeriods)=0; %#ok<*AGROW>
                    rfValsMax(e,a,channelNumber,1:numTimePeriods)=0;
                    rfValsPower(e,a,channelNumber,1:numTimePeriods)=0;
                    
                    numStimuli(e,a) = 0;
                else
                    clear erp erpBL erpST
                    erp = mean(analogData(goodPos,:),1); %#ok<*NODEF>
                    numStimuli(e,a) = length(goodPos);
                    meanLFPData(e,a,:) = erp; %#ok<*NASGU>
                    
                    for j=1:numTimePeriods
                        clear erpSegment rmsVal
                        erpSegment = erp(timePos{j});
                        rmsVal = rms(erpSegment);
                        
                        rfValsRMS(e,a,channelNumber,j) = rmsVal;
                        rfValsMax(e,a,channelNumber,j) = max(abs(erpSegment));
                        rfValsPower(e,a,channelNumber,j) = rmsVal.^2;
                    end
                end
            end
        end
        
        % Save Mean LFP data
        save(fullfile(folderOut,['meanLFPDataChan' num2str(channelNumber) fileTag '.mat']),'meanLFPData','timeVals','numStimuli');
    end
    
    % Save - numStimuli should be the same for all channels
    save(fullfile(folderOut,['rfValues' fileTag '.mat']),'rfValsRMS','rfValsMax','rfValsPower','numStimuli');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmp(measure,'Spikes')    % Case 2 - spike analysis
    
    load(fullfile(folderSegment,'LFP','lfpInfo.mat'));  % to get timeVals
    load(fullfile(folderSegment,'Spikes','spikeInfo.mat'));
    stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosGreaterThanOne);
    
    for i=1:length(neuralChannelsStored)
        channelNumber = neuralChannelsStored(i);
        SID = SourceUnitID(i);
        
        % Get Spike data
        clear spikeData neuralInfo
        load(fullfile(folderSegment,'Spikes',['elec' num2str(channelNumber) '_SID' num2str(SID) '.mat']));
        
        for a=1:aLength
            for e=1:eLength
                
                clear goodPos
                goodPos = intersect(parameterCombinations{a,e,1,end,end},stimPos);
                goodPos = setdiff(goodPos,badTrials);
                
                if isempty(goodPos)
                    rfValsRMS(e,a,channelNumber,1:numTimePeriods)=0;
                    rfValsMax(e,a,channelNumber,1:numTimePeriods)=0;
                    rfValsMean(e,a,channelNumber,1:numTimePeriods)=0;
                    
                    numStimuli(e,a) = 0;
                else
                    clear firingRate
                    [firingRate,timeValsFR] = getPSTH(spikeData(goodPos),10,timeVals);
                    numStimuli(e,a) = length(goodPos);
                    meanSpikeData(e,a,:) = firingRate;
                    
                    for j=1:numTimePeriods
                        timePos = intersect(find(timeValsFR>=timeRanges{j}(1)),find(timeValsFR<timeRanges{j}(2)));
                        clear frSegment rmsVal
                        frSegment = firingRate(timePos);
                        rmsVal = rms(frSegment);
                        
                        rfValsRMS(e,a,channelNumber,j) = rmsVal;
                        rfValsMax(e,a,channelNumber,j) = max(frSegment);
                        rfValsMean(e,a,channelNumber,j) = mean(frSegment);
                    end
                end
            end
        end
        
        % Save Mean LFP data
        save(fullfile(folderOut,['meanSpikeDataChan' num2str(channelNumber) '_SID' num2str(SID) '.mat']),'meanSpikeData','timeValsFR','numStimuli');
    end
    
    % Save - numStimuli should be the same for all channels
    save(fullfile(folderOut,'rfValues.mat'),'rfValsRMS','rfValsMax','rfValsMean','numStimuli');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmpi(measure,'CSD')    % Case 3 - CSD analysis (very similar to LFP analysis)
    
    % Get Time Ranges
    load(fullfile(folderSegment,'LFP','lfpInfo.mat'));
    timePos = cell(1,numTimePeriods);
    for i=1:numTimePeriods
        timePos{i} = intersect(find(timeVals>=timeRanges{i}(1)),find(timeVals<timeRanges{i}(2)));
    end
    
    stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosGreaterThanOne);
    
    for i=1:length(analogChannelsStored)
        channelNumber = analogChannelsStored(i);
        
        % Get CSD data
        clear signal csdData meanCSDData numStimuli
        load(fullfile(folderSegment,'CSD',['elec' num2str(channelNumber) '.mat']));
        if removeAvgRef
            csdData = csdData-avgRef;
        end
        
        for a=1:aLength
            for e=1:eLength
                
                clear goodPos
                goodPos = intersect(parameterCombinations{a,e,1,end,end},stimPos); %#ok<*USENS>
                goodPos = setdiff(goodPos,badTrials);
                
                if isempty(goodPos)
                    rfValsRMS(e,a,channelNumber,1:numTimePeriods)=0; %#ok<*AGROW>
                    rfValsMax(e,a,channelNumber,1:numTimePeriods)=0;
                    rfValsPower(e,a,channelNumber,1:numTimePeriods)=0;
                    
                    numStimuli(e,a) = 0;
                else
                    clear erp erpBL erpST
                    erp = mean(csdData(goodPos,:),1); %#ok<*NODEF>
                    numStimuli(e,a) = length(goodPos);
                    meanCSDData(e,a,:) = erp; %#ok<*NASGU>
                    
                    for j=1:numTimePeriods
                        clear erpSegment rmsVal
                        erpSegment = erp(timePos{j});
                        rmsVal = rms(erpSegment);
                        
                        rfValsRMS(e,a,channelNumber,j) = rmsVal;
                        rfValsMax(e,a,channelNumber,j) = max(abs(erpSegment));
                        rfValsPower(e,a,channelNumber,j) = rmsVal.^2;
                    end
                end
            end
        end
        
        % Save Mean LFP data
        save(fullfile(folderOut,['meanCSDDataChan' num2str(channelNumber) fileTag '.mat']),'meanCSDData','timeVals','numStimuli');
    end
    
    % Save - numStimuli should be the same for all channels
    save(fullfile(folderOut,['rfValues' fileTag '.mat']),'rfValsRMS','rfValsMax','rfValsPower','numStimuli');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmpi(measure,'Energy')    % Case 4 - Energy analysis
    
    folderMP = fullfile(folderName,'mpAnalysis');
    
    % Get Time Ranges
    clear('energyValues','downsampledFreqVals','downsampledTimeVals','numStimuli');
    load(fullfile(folderMP,['elec' num2str(goodElectrodes(1))],'energyMatlab',['mEnergy_a' num2str(1) 'e' num2str(1) '.mat']));
    
    timePos = cell(1,numTimePeriods);
    for i=1:numTimePeriods
        timePos{i} = intersect(find(downsampledTimeVals>=timeRanges{i}(1)),find(downsampledTimeVals<timeRanges{i}(2)));
    end
    numFreqPos = length(downsampledFreqVals);
    
    %stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosGreaterThanOne);

    for i=1:length(goodElectrodes)
        channelNumber = goodElectrodes(i);
        
        clear('rfValsRMS','rfValsMax','rfValsPower','numStimuliAll');
        
        for a=1:aLength
            for e=1:eLength
                
                % Get Energy data
                clear('energyValues','downsampledFreqVals','downsampledTimeVals','numStimuli');
                load(fullfile(folderMP,['elec' num2str(channelNumber)],'energyMatlab',['mEnergy_a' num2str(a) 'e' num2str(e) '.mat']));

                %clear goodPos
                %goodPos = intersect(parameterCombinations{a,e,1,end,end},stimPos); %#ok<*USENS>
                %goodPos = setdiff(goodPos,badTrials);
                
                if numStimuli==0
                    rfValsRMS(e,a,1:numTimePeriods,1:numFreqPos)=0;
                    rfValsMax(e,a,1:numTimePeriods,1:numFreqPos)=0;
                    rfValsPower(e,a,1:numTimePeriods,1:numFreqPos)=0;
                    
                    numStimuliAll(e,a) = 0;
                else
                    numStimuliAll(e,a) = numStimuli;
                    
                    for j=1:numTimePeriods
                        clear meanEnergy
                        meanEnergy = mean(energyValues(:,timePos{j}),2);

                        rfValsRMS(e,a,j,:) = sqrt(meanEnergy);
                        rfValsMax(e,a,j,:) = sqrt(max(energyValues(:,timePos{j}),[],2)); %sqrt of max energy
                        rfValsPower(e,a,j,:) = meanEnergy;
                    end
                end
            end
        end
        
        % Save Mean LFP data
        clear numStimuli
        numStimuli=numStimuliAll;
        save(fullfile(folderOut,['rfValues' num2str(channelNumber) fileTag '.mat']),'rfValsRMS','rfValsMax','rfValsPower','numStimuli','downsampledFreqVals');
    end
end
end
function stimPos = getGoodPos(subjectName,expDate,protocolName,folderSourceString,gridType,stimPosOption)

folderExtract = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData');
load(fullfile(folderExtract,'goodStimNums.mat'));
load(fullfile(folderExtract,'stimResults.mat'));

goodStimPos = stimResults.stimPosition(goodStimNums);

if exist(fullfile(folderExtract,'validStimAfterTarget.mat'),'file')
    load(fullfile(folderExtract,'validStimAfterTarget.mat'));
    if ~isempty(validStimuliAfterTarget)
        disp(['Removing ' num2str(length(validStimuliAfterTarget)) ' stimuli after target']);
    end
    goodStimPos(validStimuliAfterTarget)=-1;  % These will be not be included in either stimPos==1 or stimPos>1
end

if stimPosOption==1
    stimPos = find(goodStimPos>1);
elseif stimPosOption==2
    stimPos = find(goodStimPos==1);
else
    stimPos = find(goodStimPos>0);
end
end