% getSTA(spikeData,analogData,staTimeLims,timeVals,staLen,removeMeanSTA)
% spikeData - cell array containing spike times
% analogData - 2D array of size numTrials X stimLength
% staTimeLims - a cell array of timeLims (timesLims = [Tmin Tmax], duration
% from which spikes are taken.
% timeVals - an array of length stimLength, containing the time values
% staLen - spike triggered average duration, eg: [-0.05 0.05].
% removeMeanSTA - set to 1 if you want to remove the mean STA

function [staVals,numberOfSpikes,xsSTA] = getSTA(spikeData,analogData,staTimeLims,timeVals,staLen,removeMeanSTA)

if ~iscell(staTimeLims)
    staTimeLimsTMP = staTimeLims;
    clear staTimeLims;
    staTimeLims{1} = staTimeLimsTMP;
end

if removeMeanSTA
    analogData = analogData-repmat(mean(analogData),size(analogData,1),1);
end

Fs = round(1/(timeVals(2)-timeVals(1)));
for i=1:length(staTimeLims)
    % check if the time limit is within bounds
    timeLims = staTimeLims{i};
    if (timeLims(1)+staLen(1) < timeVals(1)) || (timeLims(2)+staLen(2) > timeVals(end))
        error('Time limit out of range');
    end 
    [staVals{i},numberOfSpikes(i)] = getSTAsingleTimeLim(spikeData,analogData,timeLims,timeVals,staLen,Fs);
end
xsSTA = staLen(1):1/Fs:staLen(2);
end

function [staVals,numberOfSpikes] = getSTAsingleTimeLim(spikeData,analogData,timeLims,timeVals,staLen,Fs)

numStim = size(analogData,1);

lim1 = staLen(1)*Fs;
lim2 = staLen(2)*Fs;
staVals = zeros(1,lim2-lim1+1);

numberOfSpikes=0;
for i=1:numStim
    clear spk 
    spk=spikeData{i};
    
    % find number of spikes in the interval
    goodSpk = spk(intersect(find(spk>=timeLims(1)),find(spk<timeLims(2))));
    
    if ~isempty(goodSpk)
        clear signal
        signal=analogData(i,:);
        
        for j=1:length(goodSpk)
            pos = max(find(timeVals<goodSpk(j)));
            numberOfSpikes=numberOfSpikes+1;
            staVals = staVals+signal(pos+lim1:pos+lim2);
        end
    end
end

if numberOfSpikes>0
    staVals=staVals/numberOfSpikes;
else
    staVals=[];
end
end