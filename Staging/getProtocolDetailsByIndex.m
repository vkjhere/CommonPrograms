function [expDate,protocolName,protocolType,baseOri,dualOri,tfRange,hasEEG] = getProtocolDetailsByIndex(index,subjectName,gridType)
% Gets the protocol details from the index. The details are retrieved from
% the protocol details file, so the protocol does not have to be extracted.
%
% Siddhesh Salelkar     14-Nov-16

if strcmpi(subjectName,'alpa')
    [expDates,protocolNames] = eval(['allProtocols' upper(subjectName(1)) subjectName(2:end) gridType]);
else
    [expDates,protocolNames] = getAllProtocols(subjectName,gridType);
end

if index > length(expDates)
    error('Invalid index.');
end

expDate = expDates{index};
protocolName = protocolNames{index};
protocolType = [];
baseOri = [];
dualOri = [];
tfRange = [];
hasEEG = [];

[~,protocolTag,protocolSet] = eval(['allProtocolDetails' upper(subjectName(1)) subjectName(2:end) gridType 'TF']);

for i=1:length(protocolSet)
    for j=1:length(protocolTag)
        [tf,do] = find(index == protocolSet{i}.indices{j});
        if ~isempty(tf) && ~isempty(do)
            protocolType = protocolTag{j};
            baseOri = protocolSet{i}.baseOri;
            if strfind(protocolTag{j},'Dual')
                dualOri = protocolSet{i}.dualOri(do);
                if strfind(protocolTag{j},'DualTF')
                    tfRange = protocolSet{i}.temporalFreqRanges{tf};
                elseif strfind(protocolTag{j},'DualTF')
                    tfRange = protocolSet{i}.temporalFreqRanges{tf};
                    centreTF = tfRange((end+1)/2);
                    tfRange = [centreTF-2 centreTF centreTF+2];
                end
            end
            hasEEG = protocolSet{i}.hasEEG;
            break;
        end
    end
end

end
