% Written by Murty Dinavahi
% Updated

function badEyeTrials = findBadTrialsFromEyeData_v2(eyeDataDeg,eyeRangeMS,FsEye,checkPeriod,bigWindowWidth,smallFixationWindowWidth)

if ~exist('bigFixationWindowWidth','var');      bigWindowWidth = 30;           end
if ~exist('smallFixationWindowWidth','var');    smallFixationWindowWidth = 5;           end

eyeDataDegX = eyeDataDeg.eyeDataDegX;
eyeDataDegY = eyeDataDeg.eyeDataDegY;
numTrials   = size(eyeDataDegX,1);

timeValsEyePos = (eyeRangeMS(1):1000/FsEye:eyeRangeMS(2)-1000/FsEye)/1000;

if ~exist('checkPeriod','var');             checkPeriod = eyeRangeMS/1000;      end
timeValsEyeCheckPos = timeValsEyePos>=checkPeriod(1) & timeValsEyePos<=checkPeriod(2);

% Find trials lying within 30
clear xTrialsBeyondFixWindow yTrialsBeyondFixWindow xTrialsNoSignals yTrialsNoSignals badEyeTrials
xTrialsBeyondBigFixWindow = sum(abs(eyeDataDegX(:,timeValsEyeCheckPos))>(bigWindowWidth/2),2);
yTrialsBeyondBigFixWindow = sum(abs(eyeDataDegY(:,timeValsEyeCheckPos))>(bigWindowWidth/2),2);
TrialsWithinBigFixWindow = find(xTrialsBeyondBigFixWindow==0 | yTrialsBeyondBigFixWindow==0 );

% do baseline correction
eyeDataDegX = eyeDataDegX(TrialsWithinBigFixWindow,:) - repmat(mean(eyeDataDegX(TrialsWithinBigFixWindow,timeValsEyeCheckPos),2),1,size(eyeDataDegX(TrialsWithinBigFixWindow,:),2));
eyeDataDegY = eyeDataDegY(TrialsWithinBigFixWindow,:) - repmat(mean(eyeDataDegY(TrialsWithinBigFixWindow,timeValsEyeCheckPos),2),1,size(eyeDataDegY(TrialsWithinBigFixWindow,:),2));

% Find bad trials
clear xTrialsBeyondFixWindow yTrialsBeyondFixWindow xTrialsNoSignals yTrialsNoSignals badEyeTrials
xTrialsBeyondSmallFixWindow = sum(abs(eyeDataDegX(:,timeValsEyeCheckPos))>(smallFixationWindowWidth/2),2);
yTrialsBeyondSmallFixWindow = sum(abs(eyeDataDegY(:,timeValsEyeCheckPos))>(smallFixationWindowWidth/2),2);
xTrialsNoSignals = sum(abs(eyeDataDegX(:,timeValsEyeCheckPos)),2);
yTrialsNoSignals = sum(abs(eyeDataDegY(:,timeValsEyeCheckPos)),2);

badEyeTrialsOutsideSmallFixWindow = xTrialsBeyondSmallFixWindow>0 | yTrialsBeyondSmallFixWindow>0 | xTrialsNoSignals==0 | yTrialsNoSignals==0;

goodEyeTrials = setdiff(TrialsWithinBigFixWindow,badEyeTrialsOutsideSmallFixWindow);
badEyeTrials = setdiff(1:numTrials,goodEyeTrials);
if badEyeTrials==0; badEyeTrials=[]; end
end