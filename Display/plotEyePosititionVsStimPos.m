% This program is used to display the eye position as a function of
% stimulus position to look for trends in eye drifts as the trial
% progresses. This program reads the data stored in
% segmentedData/eyeData/EyeDataStimPos.mat. The program
% ssaveEyeDataStimPosGRF in this folder reads the LL file and stores the
% data in the desired format

function plotEyePosititionVsStimPos(subjectName,expDate,protocolName,folderSourceString,gridType,handleX,handleY)

if ~exist('handleX','var'),          handleX=subplot(211);              end
if ~exist('handleY','var'),          handleY=subplot(212);              end

% convert to string if not already a string
if ~iscell(expDate)
    expDates{1} = expDate;
    protocolNames{1} = protocolName;
    numDays = 1;
else
    expDates=expDate;
    protocolNames=protocolName;
    numDays = length(expDates);
end

% for initialization, take the first day's data
load(fullfile(folderSourceString,'data',subjectName,gridType,expDates{1},protocolNames{1},'segmentedData','eyeData','EyeDataStimPos.mat'));
numberOfCategories = length(xs); %#ok<*USENS>

% check how many are empty
goodStimPos=[];
for i=1:numberOfCategories
    if ~isempty(xs{i})
        goodStimPos = cat(2,goodStimPos,i);
    end
end

numGoodStimPos = length(goodStimPos);
colorNames = jet(numGoodStimPos);

%initialization

for i=1:numGoodStimPos
    allDataX{i} = eyeXAllPos{goodStimPos(i)}; %#ok<*AGROW>
    allDataY{i} = eyeYAllPos{goodStimPos(i)};
    timeVals{i} = xs{goodStimPos(i)}/1000;
end

if numDays>1
    for i=2:numDays
        load(fullfile(folderSourceString,'data',subjectName,gridType,expDates{i},protocolNames{i},'segmentedData','eyeData','EyeDataStimPos.mat'));
        
        % combine data
        for j=1:numGoodStimPos
            allDataX{j} = cat(1,allDataX{j},eyeXAllPos{goodStimPos(j)});
            allDataY{j} = cat(1,allDataY{j},eyeYAllPos{goodStimPos(j)});
        end
    end
end

hold(handleX,'on'); hold(handleY,'on');
% Plot data
for i=numGoodStimPos:-1:1
    plot(handleX,timeVals{i},mean(allDataX{i}),'color',colorNames(i,:));
    plot(handleY,timeVals{i},mean(allDataY{i}),'color',colorNames(i,:));
    legendStr{numGoodStimPos-i+1} = [num2str(goodStimPos(i)) ', n=' num2str(size(allDataX{i},1))];
end

legend(handleX,legendStr,'Location','SouthEast');
axis(handleX,'tight'); axis(handleY,'tight');

axisX = axis(handleX);
axisY = axis(handleY);

% Now plot the stimulus onsets
for i=1:numGoodStimPos
    plot(handleX,(goodStimPos(i)-1)*(durationsMS.stimDurationMS+durationsMS.interStimDurationMS)/1000+zeros(1,101),axisX(3):(axisX(4)-axisX(3))/100:axisX(4),'color',colorNames(i,:));
    plot(handleY,(goodStimPos(i)-1)*(durationsMS.stimDurationMS+durationsMS.interStimDurationMS)/1000+zeros(1,101),axisY(3):(axisY(4)-axisY(3))/100:axisY(4),'color',colorNames(i,:));
end
hold(handleX,'off'); hold(handleY,'off');

set(handleX,'YTickLabel',get(handleX,'YTick'));
%set(handleY,'YTickLabel',get(handleY,'YTick'));

box(handleX,'on'); box(handleY,'on');
xlabel(handleY,'time(s)');
ylabel(handleX,'X position (deg)'); ylabel(handleY,'Y position (deg)');
title(handleX,'Eye position vs time for all conditions');
end