% Deletes a directory if it is present
function deleteDirectory(dirName)

dirName = platformSpecificName(dirName);
if exist(dirName,'dir')         
    rmdir(dirName,'s');
    disp(['deleting ' dirName]);
end
