% Deletes a file if it is present

function deleteFile(fileName)

fileName = platformSpecificName(fileName);
if exist(fileName,'file')         
    delete(fileName);         
end