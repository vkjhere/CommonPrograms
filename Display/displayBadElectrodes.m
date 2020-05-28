% This program displays bad trials and electrodes, obtained by running findBadTrialsWithEEG

function displayBadElectrodes(subjectName,expDate,protocolName,folderSourceString,gridType,capType,badTrialNameStr)

if ~exist('gridType','var');        gridType = 'EEG';                   end
if ~exist('capType','var');         capType = 'actiCap64';              end
if ~exist('badTrialNameStr','var'); badTrialNameStr = '_v5';            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderSegment = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData');
badTrialsInfo = load(fullfile(folderSegment,['badTrials' badTrialNameStr '.mat']));

%%%%%%%%%%%%%%%%%%%%%% Compare with Montage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = load([capType 'Labels.mat']); montageLabels = x.montageLabels(:,2);
x = load([capType '.mat']); montageChanlocs = x.chanlocs;

if ~isfield(badTrialsInfo,'eegElectrodeLabels') % Check if the labels match with the save labels, if these labels have been saved
    disp('Electrode labels not specified in badTrials file. Taking from Montage...');
else
    if ~isequal(montageLabels(:),badTrialsInfo.eegElectrodeLabels(:))
        error('Montage labels do not match with channel labels in badTrials');
    else
        disp('Montage labels match with channel labels in badTrials');
    end
end

if isfield(badTrialsInfo,'highPriorityElectrodeList')
    highPriorityElectrodeList = badTrialsInfo.highPriorityElectrodeList;
else
    disp('highPriorityElectrodeList not available in badTrials. Taking list from montage');
    highPriorityElectrodeList = getHighPriorityElectrodes(capType);
end
highPriorityElectrodeColor = 'g';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Topoplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot('Position',[0 0.55 0.25 0.43]);
electrodeSize = 25;
numElectrodes = length(montageLabels);

for i=1:numElectrodes
    label = num2str(i); %[num2str(i) '-' montageLabels{i}];
    montageChanlocs(i).labels = label;
end

badImpedanceElectrodes = badTrialsInfo.badElecs.badImpedanceElecs; badImpedanceElectrodeColor = 'r';
noisyElectrodes = badTrialsInfo.badElecs.noisyElecs; noisyElectrodeColor = 'm';
flatPSDElectrodes = badTrialsInfo.badElecs.flatPSDElecs; flatPSDElectrodeColor = 'b';

topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{badImpedanceElectrodes,'o',badImpedanceElectrodeColor,electrodeSize});
topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{noisyElectrodes,'o',noisyElectrodeColor,electrodeSize});
topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{flatPSDElectrodes,'o',flatPSDElectrodeColor,electrodeSize});
topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{highPriorityElectrodeList,'o',highPriorityElectrodeColor,electrodeSize});
topoplot([],montageChanlocs,'electrodes','labels','style','blank');

%%%%%%%%%%%%%%%%%%%%%%%%%% Summary Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hStats = subplot('position',[0.25 0.55 0.2 0.43]); set(hStats,'visible','off');
text(0.05,0.9,'bad Impedance Electrodes','color',badImpedanceElectrodeColor);
text(0.05,0.8,num2str(badImpedanceElectrodes(:)'),'color',badImpedanceElectrodeColor);

text(0.05,0.7,'Noisy Electrodes','color',noisyElectrodeColor);
text(0.05,0.6,num2str(noisyElectrodes(:)'),'color',noisyElectrodeColor);

text(0.05,0.5,'Flat PSD Electrodes','color',flatPSDElectrodeColor);
text(0.05,0.4,num2str(flatPSDElectrodes(:)'),'color',flatPSDElectrodeColor);

text(0.05,0.3,'High Priority Electrodes','color',highPriorityElectrodeColor);
text(0.05,0.2,num2str(highPriorityElectrodeList(:)'),'color',highPriorityElectrodeColor);

%%%%%%%%%%%%%%%%%%%%%%%% All Electrodes and Trials %%%%%%%%%%%%%%%%%%%%%%%%
badEyeTrials = badTrialsInfo.badTrialsUnique.badEyeTrials; badEyeTrialsColor = 'c';
badTrials = union(badTrialsInfo.badTrials,badEyeTrials); badTrialsColor = 'k';
numTrials = badTrialsInfo.totalTrials;
allBadTrials = badTrialsInfo.allBadTrials;
allBadTrialsMatrix = zeros(length(allBadTrials),numTrials);
for i=1:length(allBadTrials)
    if ~isnan(allBadTrials{i})
        allBadTrialsMatrix(i,allBadTrials{i}) = 1;
    end
end

h1 = getPlotHandles(1,1,[0.05 0.07 0.3 0.35]);
subplot(h1);
imagesc(1:numTrials,length(allBadTrials):-1:1,flipud(allBadTrialsMatrix),'parent',h1);
set(gca,'YDir','normal'); colormap(gray);
xlabel('# trial num','fontsize',15,'fontweight','bold');
ylabel('# electrode num','fontsize',15,'fontweight','bold');

h2 = getPlotHandles(1,1,[0.05 0.43 0.3 0.1]);
subplot(h2); cla; set(h2,'nextplot','add','XTickLabel',[]);
allBadElectrodesCount = sum(allBadTrialsMatrix,1);
stem(h2,1:numTrials,allBadElectrodesCount,'color',[0.5 0.5 0.5]); axis('tight');
ylabel('#count');

stem(h2,badTrials,sum(allBadTrialsMatrix(:,badTrials),1),'color',badTrialsColor);
stem(h2,badEyeTrials,sum(allBadTrialsMatrix(:,badEyeTrials),1),'color',badEyeTrialsColor);

%%%%%%%%%%%%%%%%%%%%%%%%%% Show electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h3 = getPlotHandles(1,1,[0.36 0.07 0.1 0.35]);
subplot(h3); cla; set(h3,'nextplot','add','XTickLabel',[]);
stem(h3,1:length(allBadTrials),sum(allBadTrialsMatrix,2),'color',[0.5 0.5 0.5]); axis('tight'); ylabel('#count');

stem(h3,badImpedanceElectrodes,sum(allBadTrialsMatrix(badImpedanceElectrodes,:),2),'color',badImpedanceElectrodeColor);
stem(h3,noisyElectrodes,sum(allBadTrialsMatrix(noisyElectrodes,:),2),'color',noisyElectrodeColor);
stem(h3,flatPSDElectrodes,sum(allBadTrialsMatrix(flatPSDElectrodes,:),2),'color',flatPSDElectrodeColor);
stem(h3,highPriorityElectrodeList,sum(allBadTrialsMatrix(highPriorityElectrodeList,:),2),'color',highPriorityElectrodeColor);
view([90 -90]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Legends %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hLegends = subplot('position',[0.36 0.43 0.1 0.1]); set(hLegends,'visible','off');

text(0.05,0.9,'badTrials (including badEye)','color',badTrialsColor);
text(0.05,0.75,'badEyeTrials','color',badEyeTrialsColor);

text(0.05,0.5,'badImpedanceElectrodes','color',badImpedanceElectrodeColor);
text(0.05,0.35,'NoisyElectrodes','color',noisyElectrodeColor);
text(0.05,0.2,'Flat PSD Electrodes','color',flatPSDElectrodeColor);
text(0.05,0.05,'HighPriorityElectrodes','color',highPriorityElectrodeColor);

%%%%%%%%%%%%% Plot raw signals for highPriorityElectrodes %%%%%%%%%%%%%%%%%
numHPElectrodes = length(highPriorityElectrodeList);
hPlots = getPlotHandles(numHPElectrodes,2,[0.5 0.05 0.45 0.9],0.05,0.01,1);
checkPeriod = [-0.500 0.750];

lfpInfo = load(fullfile(folderSegment,'LFP','lfpInfo.mat'));
timeVals = lfpInfo.timeVals;
checkPeriodIndices = timeVals>=checkPeriod(1) & timeVals<checkPeriod(2);

Fs = 1/(timeVals(2) - timeVals(1)); %Hz
params.tapers   = [3 5];
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 200];
params.trialave = 0;

for i=1:numHPElectrodes
    disp([num2str(i) ' of ' num2str(numHPElectrodes)]);
    iElec = highPriorityElectrodeList(i);
    clear x; x = load(fullfile(folderSegment,'LFP',['elec' num2str(iElec) '.mat'])); % Load EEG Data
    
    % Time domain
    plot(hPlots(i,1),timeVals,x.analogData,'g');
    hold(hPlots(i,1),'on');
    plot(hPlots(i,1),timeVals,x.analogData(badTrials,:),'r');
    xlim(hPlots(i,1),checkPeriod);
    ylabel(hPlots(i,1),[num2str(iElec) '-' montageLabels{iElec}]);
    if i<numHPElectrodes
        set(hPlots(i,1),'XTickLabel',[]);
    else
        xlabel(hPlots(i,1),'Time (seconds)');
    end
    
    % Frequency Domain
    analogDataSegment = x.analogData(:,checkPeriodIndices);
    [powerVsFreq,freqVals] = mtspectrumc(analogDataSegment',params);
    
    plot(hPlots(i,2),freqVals,log10(powerVsFreq),'g');
    hold(hPlots(i,2),'on');
    plot(hPlots(i,2),freqVals,log10(powerVsFreq(:,badTrials)),'r');
    if i<numHPElectrodes
        set(hPlots(i,2),'XTickLabel',[]);
    else
        xlabel(hPlots(i,2),'Frequency (Hz)');
    end
end

title(hPlots(1,1),['BadTrials: ' num2str(length(badTrials)) ' out of ' num2str(badTrialsInfo.totalTrials)]);
disp(['Bad Trials (n=' num2str(length(badTrials))  '): ' num2str(badTrials(:)')]);
if ~isempty(badEyeTrials)
    disp(['Bad Eye Trials (n=' num2str(length(badEyeTrials))  '): ' num2str(badEyeTrials(:)')]);
end
end
