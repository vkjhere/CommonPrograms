function closeProtocolList(dbID,dbFileName)

% if lock file does not exist, refuse to close
[dbFilePathStr,dbFileNameStr,~]=fileparts(dbFileName);
lockFile = fullfile(dbFilePathStr,['.' dbFileNameStr '.lock']);
if ~exist(lockFile,'file')
    throw(MException([mfilename ':fileNotLocked'], ...
        'File not locked, cannot close'));
end

% otherwise check if lock file matches our hostname before continuing
try
    lockFileId = fopen(lockFile,'r');
    lockingHost = fread(lockFileId);
    fclose(lockFileId);
    [cmdFail,unlockingHost] = system('hostname');
    if cmdFail
        throw(MException('system:cmdError','hostname() failed'));
    end
    lockingHost = char(lockingHost');
    if ~strcmpi(lockingHost,unlockingHost)
        throw(MException([mfilename ':hostError'], ...
            [unlockingHost ' cannot unlock file locked by ' lockingHost]));
    end
catch lockFileException
    throw(MException([mfilename ':lockFile'], ...
        ['Error verifying file lock: ' lockFileException.message]));
end

% close protocol list database file
try
    mksqlite(dbID,'close');
catch sqliteException
    rethrow(sqliteException);
end

% release lock
try
    delete(lockFile);
catch lockFileException
    throw(MException([mfilename ':lockFile'], ...
        ['Error releasing file lock: ' lockFileException.message]));
end

end
