function dataOut = combineDataGRF(dataIn)

paramsToCheck = [{'timeVals'} {'freqBL'} {'freqST'} {'freqTF'} {'timeTF'} {'frTimeVals'}]; % Need to check that these are the same for all conditions
dataOut = dataIn{1}; 

for i=1:length(paramsToCheck)
    checkDataParam(dataIn,dataOut,paramsToCheck{i});
end

% Combine data for these parameters
paramsToCombine = [{'erp'} {'SBL'} {'SST'} {'STF'} {'frVals'} {'deltaTF'} {'deltaPSD'}];
for i=1:length(paramsToCombine)
    dataOut = combineDataParam(dataIn,dataOut,paramsToCombine{i});
end
end

function checkDataParam(dataIn,dataOut,paramToCheck)
numConditions = length(dataIn);
numProtocols = length(dataOut);

for i=1:numProtocols
    paramVals = getfield(dataOut{i},paramToCheck); %#ok<*GFLD>
    for j=1:numConditions
        if ~isequal(paramVals,getfield(dataIn{j}{i},paramToCheck))
            error([paramVals ' parameters do not match']);
        end
    end
end
end
function dataOut = combineDataParam(dataIn,dataOut,paramToCombine)
numConditions = length(dataIn);
numProtocols = length(dataOut);

for i=1:numProtocols
    paramVals = getfield(dataOut{i},paramToCombine); %#ok<*GFLD>

    if strcmpi(paramToCombine,'STF') || strcmpi(paramToCombine,'deltaTF')
        paramMatrix = zeros(numProtocols,size(paramVals,1),size(paramVals,2));
        for j=1:numConditions
            paramMatrix(j,:,:) = getfield(dataIn{j}{i},paramToCombine);
        end
        dataOut{i} = setfield(dataOut{i},paramToCombine,squeeze(mean(paramMatrix,1))); %#ok<*SFLD>
    else
        paramMatrix = zeros(numProtocols,length(paramVals));
        for j=1:numConditions
            paramMatrix(j,:) = getfield(dataIn{j}{i},paramToCombine);
        end
        dataOut{i} = setfield(dataOut{i},paramToCombine,squeeze(mean(paramMatrix,1)));
    end
end
end