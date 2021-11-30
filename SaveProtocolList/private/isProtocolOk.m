function OK = isProtocolOk(protocolString)
OK = 0;
% length check
if isempty(protocolString); return; end
% format check
upos = strfind(protocolString,'_');
if isempty(upos) || upos == 1 || upos == length(protocolString); return; end
% check that protocol string is formed correctly
if ~prod(double(ismember(protocolString(1:upos-1),['A':'Z' 'a':'z']))) || ...
        ~prod(double(ismember(protocolString(upos+1:end),'0':'9'))); return; end
% check protocol number
protNum = str2num(protocolString(upos+1:end));
if isempty(protNum) || protNum < 1; return; end
% won't check for other things
OK = 1;
end
