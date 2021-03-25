function dbID = openProtocolList(dbFileName,mode)

if ~exist('mode','var'); mode = 'ro'; end

% if lock file exists, refuse to open
[dbFilePathStr,dbFileNameStr,~]=fileparts(dbFileName);
lockFile = fullfile(dbFilePathStr,['.' dbFileNameStr '.lock']);
if exist(lockFile,'file')
    throw(MException([mfilename ':fileLocked'],'File locked, cannot open'));
end

% otherwise create lock file
try
    lockFileId = fopen(lockFile,'w');
    [cmdFail,lockingHost] = system('hostname');
    if cmdFail
        throw(MException('system:cmdError','hostname() failed'));
    end
    fwrite(lockFileId,lockingHost);
    fclose(lockFileId);
catch lockFileException
    if exist(lockFile,'file'); delete(lockFile); end
    throw(MException([mfilename ':lockFile'], ...
        ['Error acquiring file lock: ' lockFileException.message]));
end

% open protocol list database file
try
    dbID = mksqlite('open',dbFileName,mode);
catch sqliteException
    rethrow(sqliteException);
end

end
