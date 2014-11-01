% comparisonData - A cell array of length L, with each cell being a 1D
% array of size Ni. This program is similar to compareAnddisplayMeanData,
% the only difference is that instead of taking in a 2D cell array and
% computing 'comparisonData' from X, we directly input comparisonData.

% The format is similar to the program compareAndDisplayData. This prgoram,
% however, takes the average of the data between the time timeLims.

function p=compareAndDisplayMeanData2(comparisonData,hPlot,hSig,testMethod,colorX,displayPlot)

defaultMethodForMoreThanTwoEntries   = 'anova';

if ~exist('hPlot','var');      hPlot = subplot(211);                    end
if ~exist('hSig','var');       hSig  = subplot(212);                    end
if ~exist('testMethod','var'); testMethod = 'default';                  end
if ~exist('displayPlot','var'); displayPlot=1;                          end

L=length(comparisonData);
if ~exist('colorX','var');     colorX = jet(L);                         end

if L==2 
    if strcmpi(testMethod,'default') || strncmpi(testMethod,'ttest',5)
        disp('Significance computed using t-test');
        [~,p,~,stats] = ttest2(comparisonData{1},comparisonData{2},0.05,'both','unequal');
        
        testParameter = 'tstat';
        testParameterValue = num2str(stats.tstat);
        testDOF = num2str(stats.df);
        pValue  = num2str(p);
        
    else
        disp('Use testMethod ttest for faster and more accurate results');
        [p,table] = anova1ForCellArrays(comparisonData,testMethod);
        
        testParameter = table{1,5};
        testParameterValue = num2str(table{2,5});
        testDOF = [num2str(table{2,3})  ', ' num2str(table{3,3})];
        pValue  = num2str(p);
        
    end
else
    if strcmpi(testMethod,'default')
        [p,table] = anova1ForCellArrays(comparisonData,defaultMethodForMoreThanTwoEntries);
    elseif strncmpi(testMethod,'ttest',5)
        disp(['Can not use t-test for more than two entries. Using ' defaultMethodForMoreThanTwoEntries ' instead']);
        [p,table] = anova1ForCellArrays(comparisonData,defaultMethodForMoreThanTwoEntries);
    else
        [p,table] = anova1ForCellArrays(comparisonData,testMethod);
    end
    
    testParameter = table{1,5};
    testParameterValue = num2str(table{2,5});
    
    if strcmp(testParameter,'F')
        testDOF = [num2str(table{2,3})  ', ' num2str(table{3,3})];
    else
        testDOF = num2str(table{2,3});
    end
    pValue  = num2str(p);
    
end

if displayPlot
    % display data
    barPlotForCellArrays(comparisonData,1:L,colorX,hPlot,0);
    
    % display the significance results
    allData=[];
    for i=1:L
        allData=[allData ; comparisonData{i}(:)]; %#ok<*AGROW>
    end
    
    tableEntries{1,1} = 'Mean';
    tableEntries{1,2} = num2str(mean(allData));
    
    tableEntries{2,1} = 'Std';
    tableEntries{2,2} = num2str(std(allData));
    
    tableEntries{3,1} = 'SEM';
    tableEntries{3,2} = num2str(std(allData)/sqrt(length(allData)));
    
    tableEntries{4,1} = testParameter;
    tableEntries{4,2} = testParameterValue;
    
    tableEntries{5,1} = 'DOF';
    tableEntries{5,2} = testDOF;
    
    tableEntries{6,1} = 'p';
    tableEntries{6,2} = pValue;
    
    makeTable(get(hSig,'Position'),tableEntries);
end
