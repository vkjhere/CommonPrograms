% This function is used to load spike and LFP data. 
% Also loads badtrials (if available) and parameterCombinations

function data = getSpikeLFPDataSingleChannel(subjectName,expDate,protocolName,folderSourceString,channelString,unitID,gridType,sideChoice,referenceChannelString,badTrialNameStr,useCommonBadTrialsFlag)

if ~exist('unitID','var');              unitID = 0;                     end
if ~exist('gridType','var');            gridType='Microelectrode';      end
if ~exist('sideChoice','var');          sideChoice=[];                  end
if ~exist('referenceChannelString','var'); referenceChannelString = ''; end
if ~exist('badTrialNameStr','var');     badTrialNameStr = '_v5';        end
if ~exist('useCommonBadTrialsFlag','var'); useCommonBadTrialsFlag = 1;  end

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
folderSegment = fullfile(folderName,'segmentedData');
folderExtract = fullfile(folderName,'extractedData');

folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');

if isnumeric(channelString) % In case channel number is entered
    channelString = ['elec' num2str(channelString)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get LFP data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear analogData
x=load(fullfile(folderLFP,channelString));
analogData=x.analogData;

% Change Reference
if strcmpi(referenceChannelString,'None') || isempty(referenceChannelString)
    % Do nothing
elseif strcmp(referenceChannelString,'AvgRef')
    disp('Changing to average reference');
    x = load(fullfile(folderLFP,'AvgRef.mat'));
    analogData = analogData - x.analogData;
else
    disp('Changing to bipolar reference');
    x = load(fullfile(folderLFP,referenceChannelString));
    analogData = analogData - x.analogData;
end
data.analogData = analogData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get timeVals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=load(fullfile(folderLFP,'lfpInfo.mat'));
data.timeVals=x.timeVals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Spike Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear spikeData
x=load(fullfile(folderSpikes,[channelString '_SID' num2str(unitID) '.mat']));
data.spikeData=x.spikeData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Get bad trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badTrialFile = fullfile(folderSegment,['badTrials' badTrialNameStr '.mat']);
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[]; allBadTrials=[];
else
    [badTrials,allBadTrials] = loadBadTrials(badTrialFile);    
end

if ~useCommonBadTrialsFlag
    badTrials = allBadTrials{str2double(channelString(5:end))};
end
data.badTrials = badTrials;
disp([num2str(length(badTrials)) ' bad trials']);

%%%%%%%%%%%%%%%%%%%%%%%% ParameterCombinations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.parameterCombinations = loadParameterCombinations(folderExtract,sideChoice);

end

function [badTrials,allBadTrials] = loadBadTrials(badTrialFile)
x=load(badTrialFile);
badTrials=x.badTrials;
allBadTrials=x.allBadTrials;
end
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract,sideChoice)

p = load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~isfield(p,'parameterCombinations2') % Not a plaid stimulus
    parameterCombinations=p.parameterCombinations;
    aValsUnique=p.aValsUnique;
    eValsUnique=p.eValsUnique;
    
    if ~isfield(p,'sValsUnique')
        sValsUnique = p.rValsUnique/3;
    else
        sValsUnique=p.sValsUnique;
    end
    
    fValsUnique=p.fValsUnique;
    oValsUnique=p.oValsUnique;
    
    if ~isfield(p,'cValsUnique')
        cValsUnique=100;
    else
        cValsUnique=p.cValsUnique;
    end
    
    if ~isfield(p,'tValsUnique')
        tValsUnique=0;
    else
        tValsUnique=p.tValsUnique;
    end 
else
    [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
        fValsUnique,oValsUnique,cValsUnique,tValsUnique] = makeCombinedParameterCombinations(folderExtract,sideChoice);
end
end