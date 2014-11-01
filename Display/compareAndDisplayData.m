% This program is different from the program compareAndDisplayData.m in the
% parent folder. This program is an extension that allows a comparison of
% more than two datasets. When only two datasets are presents, it uses the
% old program (which is now called compareAnddisplayData_ttest).

% X - A cell array of length L, with each cell being a 2D array of size
% NixK (i=1...L). Ni is the number of stimuli. xs must be of length K.

function p=compareAndDisplayData(X,xs,hPlot,hSig,testMethod,colorX,colorSig,displayPlot,showStdErr)

defaultMethodForMoreThanTwoEntries   = 'anova';

if ~exist('hPlot','var');      hPlot = subplot(211);                     end
if ~exist('hSig','var');       hSig  = subplot(212);                     end
if ~exist('testMethod','var'); testMethod = 'default';                   end
if ~exist('colorSig','var');   colorSig = 'k';                           end
if ~exist('displayPlot','var'); displayPlot=1;                           end
if ~exist('showStdErr','var');   showStdErr = 0;                         end

L=length(X);
if ~exist('colorX','var');     colorX = jet(L);                          end
if ~exist('xs','var');          xs = 1:size(X{1},2);                     end

if L==2 
    if strcmpi(testMethod,'default') || strncmpi(testMethod,'ttest',5)
        disp('Significance computed using t-test');
        p=compareAndDisplayData_ttest(X{1},X{2},xs,hPlot,hSig,colorX(1,:),colorX(2,:),colorSig,displayPlot,showStdErr);
    else
        disp('Use testMethod ttest for faster and more accurate results');
        p=compareAndDisplayData_AnovaOrKW(X,xs,hPlot,hSig,testMethod,colorX,colorSig,displayPlot,showStdErr);
    end
else
    if strcmpi(testMethod,'default')
        p=compareAndDisplayData_AnovaOrKW(X,xs,hPlot,hSig,defaultMethodForMoreThanTwoEntries,colorX,colorSig,displayPlot,showStdErr);
    elseif strncmpi(testMethod,'ttest',5)
        disp(['Can not use t-test for more than two entries. Using ' defaultMethodForMoreThanTwoEntries ' instead']);
        p=compareAndDisplayData_AnovaOrKW(X,xs,hPlot,hSig,defaultMethodForMoreThanTwoEntries,colorX,colorSig,displayPlot,showStdErr);
    else
        p=compareAndDisplayData_AnovaOrKW(X,xs,hPlot,hSig,testMethod,colorX,colorSig,displayPlot,showStdErr);
    end
end

axis(hPlot,'tight');
axis(hSig,'tight');