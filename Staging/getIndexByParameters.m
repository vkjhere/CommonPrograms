function indices = getIndexByParameters(subjectName,gridType,...
    azimuth,elevation,sigma,spatialFreq,orientation,contrast,temporalFreq,...
    folderSourceString,protocolPrefix)
% Gets the indices of extracted protocols matching the given parameters.
% Empty parameters are ignored for computing the match.
%
% Siddhesh Salelkar     15-Sep-16

if ~exist('folderSourceString','var'); folderSourceString = 'D:'; end
if ~exist('protocolPrefix','var'); protocolPrefix = 'GRF_'; end

if strcmpi(subjectName,'abu') || strcmpi(subjectName,'rafiki') || strcmpi(subjectName,'alpa')
    [expDates,protocolNames] = eval(['allProtocols' upper(subjectName(1)) subjectName(2:end)  gridType]);
else
    [expDates,protocolNames] = getAllProtocols(subjectName,gridType);
end

indices = [];
for i=1:length(expDates)
    if strcmp(protocolNames{i}(1:length(protocolPrefix)),protocolPrefix)
        [a,e,s,f,o,c,t,paramsFound] = loadParameterCombinations(subjectName,expDates{i},protocolNames{i},gridType,folderSourceString);
        if ~paramsFound
            continue;
        end
        clear indexMatch
        indexMatch = ones(1,7);
        if ~isempty(azimuth) && ~prod(ismember(azimuth,a))
            indexMatch(1) = 0;
        end
        if ~isempty(elevation) && ~prod(ismember(elevation,e))
            indexMatch(2) = 0;
        end
        if ~isempty(sigma) && ~prod(ismember(sigma,s))
            indexMatch(3) = 0;
        end
        if ~isempty(spatialFreq) && ~prod(ismember(spatialFreq,f))
            indexMatch(4) = 0;
        end
        if ~isempty(orientation) && ~prod(ismember(orientation,o))
            indexMatch(5) = 0;
        end
        if ~isempty(contrast) && ~prod(ismember(contrast,c))
            indexMatch(6) = 0;
        end
        if ~isempty(temporalFreq) && ~prod(ismember(temporalFreq,t))
            indexMatch(7) = 0;
        end
        if prod(indexMatch)
            indices = [indices i]; %#ok<AGROW>
        end
    end
end

end

function [aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique,paramsFound] = loadParameterCombinations(subjectName,expDate,protocolName,gridType,folderSourceString)

aValsUnique = []; eValsUnique = []; sValsUnique = []; fValsUnique = [];
oValsUnique = []; cValsUnique = []; tValsUnique = []; paramsFound = 0;

pcFileName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat');
if ~exist(pcFileName,'file')
    return;
end

load(pcFileName);
paramsFound = 1;
end
