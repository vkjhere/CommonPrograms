function insertIntoProtocolList(dbFileName,expDate,protocolName,stimType,comment)

% convert from ddmmyy to yymmdd while saving: this is required to use the
% ORDER BY clause when retrieving records
expDate = [expDate(5:6) expDate(3:4) expDate(1:2)];

try
    % open protocol list database file
    dbID = openProtocolList(dbFileName,'rw');
    % insert into protocol list table in database
    mksqlite(dbID, ...
        ['INSERT INTO PROTOCOLLIST (EXPDATE,PROTOCOLNAME,STIMTYPE,COMMENT) ' ...
        'VALUES (''' expDate ''',''' protocolName ''',' num2str(stimType) ',''' comment ''');']);
    % close database
    closeProtocolList(dbID,dbFileName);
catch sqliteException
    if exist('dbID','var'); closeProtocolList(dbID,dbFileName); end
    rethrow(sqliteException);
end

end
