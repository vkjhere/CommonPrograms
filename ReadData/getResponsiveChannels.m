function responsiveChannels = getResponsiveChannels(subjectName,expDate,protocolName,gridType,folderSourceString,...
    blPeriod,stPeriod,firingThresh,frModThresh,alpha,displayDataFlag,saveDataFlag)
% Get the responsive channels for a particular protocol, based on spiking
% activity during baseline versus stimulus. Channels are filtered based on
% a combination of absolute firing threshold (firingThresh, spk/sec),
% firing modulation threshold (frModThresh = (stFR-blFR)/(stFR+blFR)) and
% parameteric (ANOVA) or non-parametric (Sign Rank) measures of difference
% between baseline and stimulus firing activity.
%
% Siddhesh Salelkar 05-Sep-2016

if ~exist('blPeriod','var'); blPeriod = [-1.5 0]; end
if ~exist('stPeriod','var'); stPeriod = [0 1.5]; end
if ~exist('firingThresh','var'); firingThresh = 5; end % conservative!
if ~exist('frModThresh','var'); frModThresh = 0.2; end
if ~exist('alpha','var'); alpha = 0.05; end
if ~exist('displayDataFlag','var'); displayDataFlag = 1; end
if ~exist('saveDataFlag','var'); saveDataFlag = 0; end

folderData = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderExtract = fullfile(folderData,'extractedData');
folderSegment = fullfile(folderData,'segmentedData');

load(fullfile(folderExtract,'parameterCombinations.mat'));
load(fullfile(folderSegment,'Spikes','spikeInfo.mat'));

badTrialFile = fullfile(folderSegment,'badTrials.mat');
if ~exist(badTrialFile,'file')
    badTrials = [];
    disp('Bad trial file does not exist');
else
    load(badTrialFile,'badTrials');
end

%% Get list of channels responsive for each parameter combination
% Since each channel can potentially fire on only a certain combination of
% parameters, we find the list of channels separately for each parameter
% combination and then pool them.

responsiveChannelsStim = cell(length(aValsUnique),length(eValsUnique),length(sValsUnique), ...
    length(fValsUnique),length(oValsUnique),length(cValsUnique),length(tValsUnique));

for a = 1:length(aValsUnique)
for e = 1:length(eValsUnique)
for s = 1:length(sValsUnique)
for f = 1:length(fValsUnique)
for o = 1:length(oValsUnique)
for c = 1:length(cValsUnique)
for t = 1:length(tValsUnique)

    trialNumbers = parameterCombinations{a,e,s,f,o,c,t}; %#ok<*USENS>
    trialNumbers = setdiff(trialNumbers,badTrials); % could also use electrode-specfic bad trials...
    
%     disp(['a=' num2str(aValsUnique(a)) ',e=' num2str(eValsUnique(e)) ',s=' sValsUnique(s) ...
%         ',f=' num2str(fValsUnique(f)) ',o=' num2str(oValsUnique(o)) ',c=' num2str(cValsUnique(c)) ...
%         ',t=' num2str(tValsUnique(t)) ': n=' num2str(length(trialNumbers)) ' trials']);

    channelArray = [];
    chanCount = 0;
    
    for chan = 1:length(neuralChannelsStored)
        clear spikingActivity;
        spikingActivity = load(fullfile(folderSegment,'Spikes',['elec' num2str(neuralChannelsStored(chan)) '_SID0.mat']));

        % separate spikes into baseline and stimulus
        numSpikesBL = cellfun(@(x) length(find(x >= blPeriod(1) & x < blPeriod(2))), ...
            spikingActivity.spikeData(trialNumbers)) / diff(blPeriod);
        numSpikesST = cellfun(@(x) length(find(x >= stPeriod(1) & x < stPeriod(2))), ...
            spikingActivity.spikeData(trialNumbers)) / diff(stPeriod);

        % mean spiking activity between BL and ST is significantly different?
        pval = anova1([numSpikesBL(:) numSpikesST(:)],[],'off');
        
        % median change in spiking activity from BL to ST is significantly different from 0?
        pval2 = signrank(numSpikesST,numSpikesBL);

        % calculate firing rates and modulation
        blFR = mean(numSpikesBL);
        stFR = mean(numSpikesST);
        frMod = (stFR - blFR) / (stFR + blFR);
        excited = blFR < firingThresh && stFR >= firingThresh;
        inhibited = blFR >= firingThresh && stFR < firingThresh;

        % check if it meets all the criteria
        if (excited || inhibited) && abs(frMod) >= frModThresh && (pval < alpha || pval2 < alpha)
            chanCount = chanCount + 1;
            channelArray(chanCount).num = neuralChannelsStored(chan); %#ok<*AGROW>
            channelArray(chanCount).blFR = blFR;
            channelArray(chanCount).stFR = stFR;
            channelArray(chanCount).frMod = frMod;
            channelArray(chanCount).excited = double(excited);
            channelArray(chanCount).inhibited = double(inhibited);
            channelArray(chanCount).panova = pval;
            channelArray(chanCount).psignrank = pval2;
        end
    end

%     if ~isempty(channelArray); disp([channelArray.num]); end
    responsiveChannelsStim{a,e,s,f,o,c,t} = channelArray;
end
end
end
end
end
end
end

%% Get all responsive channels with their properties

allChannels = responsiveChannelsStim(:);
responsiveChanList = [];
responsiveChanFrModList = [];
blFRList = []; stFRList = [];
excitedList = []; inhibitedList = [];
panovaList = []; psignrankList = [];
for i=1:length(allChannels)
    if ~isempty(allChannels{i})
        responsiveChanList = [responsiveChanList [allChannels{i}.num]];
        responsiveChanFrModList = [responsiveChanFrModList [allChannels{i}.frMod]];
        blFRList = [blFRList [allChannels{i}.blFR]];
        stFRList = [stFRList [allChannels{i}.stFR]];
        excitedList = [excitedList [allChannels{i}.excited]];
        inhibitedList = [inhibitedList [allChannels{i}.inhibited]];
        panovaList = [panovaList [allChannels{i}.panova]];
        psignrankList = [psignrankList [allChannels{i}.psignrank]];
    end
end

% summary
responsiveChannels = unique(responsiveChanList);
disp([num2str(length(responsiveChanList)) ' electrode(s) responsive, ' num2str(length(responsiveChannels)) ' unique']);

%% Display data and save if required

if displayDataFlag
    fprintf('\nFiring rate threshold = %.2f spk/s\n',firingThresh);
    fprintf('Firing rate modulation index threshold = %.2f\n',frModThresh);
    fprintf('Alpha threshold = %.4f\n\n',alpha);
    
    fprintf('Channel   FrMod   Response     ANOVA     SignRank   MeanFR   SNR\n');
    fprintf('-------   -----   --------     -----     --------   ------   ---\n');
    
    % display in console
    for i=1:length(responsiveChannels)
        respIndices = find(responsiveChanList == responsiveChannels(i));
        frModChan(i) = mean(responsiveChanFrModList(respIndices));
        frModChanSE(i) = std(responsiveChanFrModList(respIndices))/sqrt(length(respIndices));
        for j = 1:length(respIndices)
            if excitedList(respIndices(j))
                respString = 'excited  ';
                meanFR = num2str(stFRList(respIndices(j)),'%6.2f');
            elseif inhibitedList(respIndices(j))
                respString = 'inhibited';
                meanFR = num2str(blFRList(respIndices(j)),'%6.2f');
            else
                respString = 'none     ';
                meanFR = ' ---- ';
            end
            if j == 1
                fprintf('%3d       %5.2f   %s    %6.4f    %6.4f     %s     %.2f\n', ...
                    responsiveChannels(i),responsiveChanFrModList(respIndices(j)), ...
                    respString,panovaList(respIndices(j)),psignrankList(respIndices(j)), ...
                    meanFR,channelSNR(folderSegment,responsiveChannels(i)));
            else
                fprintf('          %5.2f   %s    %6.4f    %6.4f     %s\n', ...
                    responsiveChanFrModList(respIndices(j)),respString, ...
                    panovaList(respIndices(j)),psignrankList(respIndices(j)),meanFR);
            end
        end
    end
    
    % plot bar graph
    bar(frModChan); hold on; errorbar(frModChan,frModChanSE,'LineStyle','none');
    xlim([0 length(responsiveChannels)+1]); ylim([-1 1]);
    set(gca,'XTick',0:length(responsiveChannels)+1);
    set(gca,'XTickLabel',[' ' convertArrayToCellArray(responsiveChannels) ' ']);
    xlabel('Electrodes'); ylabel('Firing rate modulation index');
end

if saveDataFlag
    save(fullfile(folderData,'responsiveChannels'), ...
        'responsiveChannelsStim','responsiveChannels', ...
        'blPeriod','stPeriod','firingThresh','frModThresh','alpha');
end

end

function snrval = channelSNR(folderSegment,channel)
    
    x = load(fullfile(folderSegment,'Segments',['elec' num2str(channel) '.mat']));
    snrval = getSNR(x.segmentData');
end
