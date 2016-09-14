function arrayOut = convertArrayToCellArray(arrayIn,formatSpec,prefix,suffix)
% Converts a numeric array into a cell array of strings
%
% Siddhesh Salelkar     14-Sep-2016

if ~exist('formatSpec','var'); formatSpec = '%d'; end
if ~exist('prefix','var'); prefix = ''; end
if ~exist('suffix','var'); suffix = ''; end

tempArray = arrayIn(:);
for i=1:length(tempArray)
    arrayOut{i} = [prefix num2str(tempArray(i),formatSpec) suffix]; %#ok<*AGROW>
end
arrayOut = reshape(arrayOut,size(arrayIn));

end
