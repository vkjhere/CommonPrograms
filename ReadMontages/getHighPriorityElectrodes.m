% high priority electrodes are electrodes which are given more importance when
% determining bad trials. Usually they are the visual electrodes.

% If electrode group is not empty, the the highPriorityList is one of the
% electrodeGroups specified in electrodePositionOnGrid

function highPriorityElectrodeNums = getHighPriorityElectrodes(capType,electrodeGroup)
if ~exist('electrodeGroup','var');  electrodeGroup='';                  end

[~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(1,'EEG',[],capType);

if ~isempty(electrodeGroup)
    highPriorityElectrodeNums = cell2mat(electrodeGroupList(strcmp(groupNameList,electrodeGroup)));
end
end