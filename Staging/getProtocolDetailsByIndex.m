function [expDates,protocolNames,protocolTypes,baseOris,dualOris,tfRanges] = getProtocolDetailsByIndex(index,subjectName,gridType)
% Gets the protocol details from the index. The details are retrieved from
% the protocol details file, so the protocol does not have to be extracted.
%
% Siddhesh Salelkar     14-Nov-16

if strcmpi(subjectName,'alpa')
    [allExpDates,allProtocolNames] = eval(['allProtocols' upper(subjectName(1)) subjectName(2:end) gridType]);
else
    [allExpDates,allProtocolNames] = getAllProtocols(subjectName,gridType);
end

if index > length(allExpDates)
    error('Invalid index.');
end

expDates = cell(1,length(index));
protocolNames = cell(1,length(index));
protocolTypes = cell(1,length(index));
baseOris = cell(1,length(index));
dualOris = cell(1,length(index));
tfRanges = cell(1,length(index));

[~,protocolTag,protocolSet] = eval(['allProtocolDetails' upper(subjectName(1)) subjectName(2:end) gridType 'TF']);

for ii=1:length(index)
    expDates{ii} = allExpDates{index(ii)};
    protocolNames{ii} = allProtocolNames{index(ii)};
    if nargout == 2; continue; end

    for i=1:length(protocolSet)
        for j=1:length(protocolTag)
            [tf,do] = find(index(ii) == protocolSet{i}.indices{j});
            if ~isempty(tf) && ~isempty(do)
                protocolTypes{ii} = protocolTag{j};
                baseOris{ii} = protocolSet{i}.baseOri;
                if strfind(protocolTag{j},'Dual')
                    dualOris{ii} = protocolSet{i}.dualOri(do);
                    if strfind(protocolTag{j},'DualTF')
                        tfRanges{ii} = protocolSet{i}.temporalFreqRanges{tf};
                    elseif strfind(protocolTag{j},'DualCon')
                        tfRange = protocolSet{i}.temporalFreqRanges{tf};
                        centreTF = tfRange((end+1)/2);
                        tfRanges{ii} = [centreTF-2 centreTF centreTF+2];
                    end
                end
                break;
            end
        end
    end
end

end
