function highPriorityElectrodeNums = getHighPriorityElectrodes(capType)
% high priority electrodes are electrodes which are given more importance when
% determining bad trials. Usually they are the visual electrodes.

if strcmp(capType,'actiCap64')
    highPriorityElectrodeNums = [24 26 29 30 31 57 58 61 62 63];
elseif strcmp(capType,'actiCap31Posterior')
    highPriorityElectrodeNums = [15 16 18 19 23 24 25 28 29 30];
else
    error('capType not found');
end
end