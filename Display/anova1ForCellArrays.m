function [p,table] = anova1ForCellArrays(X,testMethod)

comparisonData = []; groupValue = [];
L = length(X);

for j=1:L
    comparisonData = [comparisonData ; X{j}(:)]; %#ok<*AGROW>
    groupValue     = [groupValue ; j+zeros(length(X{j}),1)];
end

if strcmpi(testMethod,'anova')
    [p,table] = anova1(comparisonData,groupValue,'off');
elseif strcmpi(testMethod,'kruskalWallis') || strcmpi(testMethod,'KW')
    [p,table] = kruskalwallis(comparisonData,groupValue,'off');
end
    
