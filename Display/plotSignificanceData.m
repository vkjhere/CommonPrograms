% This program plots the p-values, two lines at pCutOff and 0.05/length(p)
% and shows the values for which the p values reach significance

function plotSignificanceData(p,xs,hSig,colorSig,pCutOff)

if ~exist('hSig','var');       hSig  = subplot(212);                    end
if ~exist('colorSig','var');   colorSig = 'k';                          end
if ~exist('pCutOff','var');    pCutOff=0.01;                            end

color1 = 'g'; % 0.01
color2  = 'r'; % Bonferroni

K=length(p);
% Plot data
if ~exist('xs','var');          xs = 1:K;                               end

%Plot Significance
set(hSig,'Nextplot','replace');
plot(hSig,xs,log10(p),'Color',colorSig);
set(hSig,'Nextplot','add');
plot(hSig,xs,log10(pCutOff+zeros(1,K)),'Color',color1);
plot(hSig,xs,log10(0.05/K+zeros(1,K)),'Color',color2);

% Significant points
axis(hSig,'tight');
axisLims = axis(hSig);
axis(hSig,[axisLims(1:2) axisLims(3) 0]);

sig1 = find(p<pCutOff);
sig2 = find(p<0.05/K);

if ~isempty(sig1)
    plot(hSig,xs(sig1),0,'Color',color1,'marker','*');
end

if ~isempty(sig2)
    plot(hSig,xs(sig2),0,'Color',color2,'marker','*');
end
set(hSig,'Nextplot','replace');

end