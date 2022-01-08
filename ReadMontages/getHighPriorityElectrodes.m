function highPriorityElectrodeNums = getHighPriorityElectrodes(capType,gridType,electrodeGroup)
% high priority electrodes are electrodes which are given more importance when
% determining bad trials. Usually they are the visual electrodes.
if strcmp(capType,'actiCap64')
    [~,~,~,electrodeGroupList,groupNameList] = electrodePositionOnGrid(1,gridType,[],2);
    highPriorityElectrodeNums = cell2mat(electrodeGroupList(strcmp(groupNameList,electrodeGroup)));
elseif strcmp(capType,'actiCap31Posterior')
    highPriorityElectrodeNums = [15 16 18 19 23 24 25 28 29 30];
elseif strcmp(capType,'actiCap64_2019')
    [~,~,~,electrodeGroupList,groupNameList] = electrodePositionOnGrid(1,gridType,[],6);
    highPriorityElectrodeNums = cell2mat(electrodeGroupList(strcmp(groupNameList,electrodeGroup)));
else
    error('capType not found');
end
end