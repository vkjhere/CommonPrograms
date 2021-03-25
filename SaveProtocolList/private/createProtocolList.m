function createProtocolList(dbFileName)

% check if database file already exists
if exist(dbFileName,'file')
    throw(MException([mfilename ':fileExists'], ...
        'Cannot create database, file already exists'));
end

try
    % create protocol list database file
    dbID = openProtocolList(dbFileName,'rwc');
    % initialize protocol list table in database
    mksqlite(dbID, ...
        ['CREATE TABLE PROTOCOLLIST(' ...
            'EXPDATE VARCHAR NOT NULL, ' ...
            'PROTOCOLNAME VARCHAR NOT NULL, ', ...
            'STIMTYPE INT NOT NULL, ' ...
            'COMMENT VARCHAR, ' ...
            'PRIMARY KEY(EXPDATE,PROTOCOLNAME)' ...
        ');']);
    % close database
    closeProtocolList(dbID,dbFileName);
catch sqliteException
    if exist('dbID','var'); closeProtocolList(dbID,dbFileName); end
    rethrow(sqliteException);
end

end
