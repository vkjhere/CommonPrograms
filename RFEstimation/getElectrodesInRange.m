function electrodeNums = getElectrodesInRange(subjectName,expDate,protocolName,folderSourceString,dRange,stimulusCenter,useTheseElectrodes,gridType)

if ~exist('stimulusCenter','var');       stimulusCenter=[];             end
if ~exist('useTheseElectrodes','var');   useTheseElectrodes=[];         end
if ~exist('gridType','var');             gridType = 'Microelectrode';   end

folderRFData = fullfile(removeIfPresent(fileparts(mfilename('fullpath')),fullfile('ProgramsMAP','CommonPrograms','RFEstimation')),'DataMAP','ReceptiveFieldData');

% get Stimulus position information
if isempty(stimulusCenter)
    load(fullfile(folderSourceString,subjectName,expDate,protocolName,extractedData,'parameterCombinations.mat'));
else
    aValsUnique=stimulusCenter(1);
    eValsUnique=stimulusCenter(2);
end

% get RF information
load(fullfile(folderRFData,[subjectName gridType 'RFData.mat']));

count=1;
if isempty(useTheseElectrodes)
    useTheseElectrodes=highRMSElectrodes;
end

for i=1:length(useTheseElectrodes)
    azi = rfStats(useTheseElectrodes(i)).meanAzi;
    ele = rfStats(useTheseElectrodes(i)).meanEle;
    
    d = sqrt(sum((azi-aValsUnique)^2+(ele-eValsUnique)^2));
    
    if d>=dRange(1) && d<dRange(2)
        electrodeNums(count) = useTheseElectrodes(i); %#ok<AGROW>
        count=count+1;
    end
end

if count==1 % Nothing found
    electrodeNums=[];
end
end