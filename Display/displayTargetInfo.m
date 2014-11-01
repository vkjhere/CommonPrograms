function displayTargetInfo(gridPositionPerformanceVsTargetPos,gridPositionPercentTargets,gridPositionCorrect,monkeyName,expDates,protocolNames,folderSourceString,gridType,removeBreaks,removeIgnores)

if iscell(expDates)
    allTargetPos=[];
    allCatchTrials = [];
    allEOTCodes = [];
    
    goodTargetPos = [];
    goodStimPos   = [];
    
    for i=1:length(expDates)
        clear expDate protocolName
        expDate=expDates{i};
        protocolName=protocolNames{i};
        
        clear allTrials goodTrials stimData 
        load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'extractedData','BehaviorData.mat'));
        
        % concatenate
        allTargetPos = cat(2,allTargetPos,allTrials.targetPosAllTrials);
        allCatchTrials = cat(2,allCatchTrials,allTrials.catchTrials);
        allEOTCodes  = cat(2,allEOTCodes,allTrials.eotCodes);
        
        goodTargetPos = cat(2,goodTargetPos,goodTrials.targetPos);
        goodStimPos   = cat(2,goodStimPos,stimData.stimPos);
    end
else
    expDate=expDates;
    protocolName=protocolNames;
    load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'extractedData','BehaviorData.mat')); % Single day
    
    allTargetPos = [allTrials.targetPosAllTrials];
    allCatchTrials = [allTrials.catchTrials];
    allEOTCodes  = [allTrials.eotCodes];
    goodTargetPos = goodTrials.targetPos;
    goodStimPos   = stimData.stimPos;
end

% Assign a target position of maxTargetPos+1 to all catch trials
if sum(allCatchTrials) == length(allCatchTrials) % All Trials are catch trials
    maxTargetPos = max(allTargetPos);
else
    maxTargetPos = max(allTargetPos(allCatchTrials==0));
end
allTargetPos(allCatchTrials==1)=maxTargetPos+1;

gridWidthLL = gridPositionPerformanceVsTargetPos(3)/(maxTargetPos+1);
gridHeightLL = gridPositionPerformanceVsTargetPos(4);

targetsAtPos = zeros(1,maxTargetPos+1);
for i=1:maxTargetPos+1
    tmpEOTCodes = allEOTCodes(allTargetPos==i);
    %targetsAtPos(i) = length(tmpEOTCodes);
    gridOrigin = [gridPositionPerformanceVsTargetPos(1)+ (i-1)*gridWidthLL gridPositionPerformanceVsTargetPos(2)];
    gridPositionLLBar = [gridOrigin gridWidthLL gridHeightLL];

    [~,~,targetsAtPos(i)] = displayEOTCodes(gridPositionLLBar,tmpEOTCodes,protocolName,1,[],removeBreaks,removeIgnores);
end

subplot('Position',gridPositionPercentTargets);
targetAtPosPercent = 100*targetsAtPos/sum(targetsAtPos);
stem(targetAtPosPercent); 
axis([0.5 maxTargetPos+1.5 0 max(targetAtPosPercent)]);
title(['N=' num2str(sum(targetsAtPos)) ' (rightmost=catch)']);
ylabel('Percent');

%%%%%%%%%%%%%%%%%%%%%%%%%%% Target Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot('Position',gridPositionCorrect{1});
numGoodTrials = length(goodTargetPos);

goodTargetPosVal = zeros(1,max(goodTargetPos));
for i=1:max(goodTargetPos)
    goodTargetPosVal(i) = 100*length(find(goodTargetPos==i))/numGoodTrials;
end
stem(goodTargetPosVal);
title(['Target Pos of N=' num2str(numGoodTrials) ' good trials']);

%%%%%%%%%%%%%%%%%%%%%%%%% Stimulus Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot('Position',gridPositionCorrect{2});
numGoodStims = length(goodStimPos);

goodStimPosVal = zeros(1,max(goodStimPos));
for i=1:max(goodStimPos)
    goodStimPosVal(i) = 100*length(find(goodStimPos==i))/numGoodStims;
end
stem(goodStimPosVal);
title(['Pos of N=' num2str(numGoodStims) ' good stimuli']);
end
