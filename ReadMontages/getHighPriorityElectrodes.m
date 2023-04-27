% high priority electrodes are electrodes which are given more importance when
% determining bad trials. Usually they are the visual electrodes.

% If electrode group is not empty, the the highPriorityList can be as many
% of the specified electrodeGroups  

function highPriorityElectrodeNums = getHighPriorityElectrodes(capType,electrodeGroup)
if ~exist('electrodeGroup','var');  electrodeGroup='';                  end
[~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(1,'EEG',[],capType);

if ~isempty(electrodeGroup)
    if iscell(electrodeGroup)
        highPriorityElectrodeNums = [];
        for ilen = 1:length(electrodeGroup)
            highPriorityElectrodeNumsEachGroup = cell2mat(electrodeGroupList(strcmp(groupNameList,electrodeGroup{ilen})));
            highPriorityElectrodeNums = [highPriorityElectrodeNums,highPriorityElectrodeNumsEachGroup];
        end
    else
        highPriorityElectrodeNums = cell2mat(electrodeGroupList(strcmp(groupNameList,electrodeGroup)));
    end
end
end