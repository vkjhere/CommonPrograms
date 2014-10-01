function [Y,xs] = convertSpikeTimes2Bins(X,tRangeS,dMS)

if ~exist('dMS','var');                    dMS=1;                       end

dS = dMS/1000;
xs = tRangeS(1)+dS/2:dS:tRangeS(2);
numBins = length(xs);
numStimuli = length(X);

Y = zeros(numBins,numStimuli);

for i=1:numStimuli
    tVals=X{i};
    y = tVals(intersect(find(tVals>=tRangeS(1)),find(tVals<tRangeS(2))))-tRangeS(1); % spike times normalized to 0
    Y(max(1,min(unique(ceil(y/dS)),numBins)),i)=1;
end