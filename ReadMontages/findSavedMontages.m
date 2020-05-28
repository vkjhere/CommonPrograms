% This program shows all saved Montages, which are kept in the folder
% 'Montages' that must be in Matlab path. Further, the actiCap64 layout
% file actiCap64.mat must be present under the folder
% Montages\Layouts\actiCap64. This program finds the location of the
% Layouts folder based on this information and then other layouts that are
% available in that folder.

function montageFolderNames = findSavedMontages

layoutsFolder = fileparts(fileparts(which('actiCap64.mat')));

if isempty(layoutsFolder)
    error('actiCap64.mat file not found in the search path');
else
    allFiles = dir(layoutsFolder);
    
    count = 1;
    for i=1:length(allFiles)
        thisFile = allFiles(i);
        if (thisFile.isdir) && ~strcmp(thisFile.name,'.') && ~strcmp(thisFile.name,'..')
            montageFolderNames{count} = thisFile.name; %#ok<AGROW>
            count=count+1;
        end
    end
end