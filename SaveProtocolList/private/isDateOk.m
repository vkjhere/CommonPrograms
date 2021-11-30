function OK = isDateOk(dateString)
OK = 0;
% length check
if length(dateString) ~= 6; return; end
% check fields are all numbers
dd = str2num(dateString(1:2));
mm = str2num(dateString(3:4));
yy = str2num(dateString(5:6));
if isempty(dd) || isempty(mm) || isempty(yy); return; end
% check dd and mm are within permissible limits
if yy < 1 || dd < 1 || dd > 31 || mm < 1 || mm > 12; return; end
% won't check for validity of 28/29/30/31 day months
OK = 1;
end
