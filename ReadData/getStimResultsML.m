% This program generates a dummy stimResults file. In monkeyLogic (ML),
% each stimulus is identified by a particular number. These numbers are
% simply mapped to the orientation, while the rest of the parameters are
% simply set to zero.

function getStimResultsML(folderExtract,stimNumbers)

dummyList = zeros(1,length(stimNumbers));
stimResults.azimuth = dummyList;
stimResults.elevation = dummyList;
stimResults.contrast = dummyList;
stimResults.temporalFrequency = dummyList;
stimResults.radius = dummyList;
stimResults.sigma = dummyList;
stimResults.orientation = stimNumbers;
stimResults.spatialFrequency = dummyList;
stimResults.side = 0;

save(fullfile(folderExtract,'stimResults.mat'),'stimResults');
end