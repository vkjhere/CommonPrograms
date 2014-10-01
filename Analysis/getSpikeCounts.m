% H =getSpikeCounts(X,tRange)
% This function returns the number of spikes in each trial of the cell
% array X between times tRange (in seconds)

% Supratim Ray 09/29/09
% Modified: 01/03/11. tRange is now a vector

function H = getSpikeCounts(X,tRange)

if isempty(X)
    H=0;
else
    numTrials = length(X);
    
    for i=1:numTrials
        spk = X{i};
        H(i) = length(intersect(find(spk>=tRange(1)),find(spk<tRange(2))));
    end
end