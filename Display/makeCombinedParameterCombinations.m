% When both Gabor0 and Gabor1 stimuli are used together, say in case of a
% plaid, it is possible to combine the parameterCombinations when only one
% out of two gabor paramaters takes more than 1 value. This program will not
% work if both parameters have more than one value. For example, if grating
% 1 has 4 orientations and 1 TF while grating2 has 1 orientation and 10
% TFs, this program will return a combined parameterCombinations matrix
% with 4 orientations and 10 TFs.

function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = makeCombinedParameterCombinations(folderExtract)

p = load(fullfile(folderExtract,'parameterCombinations.mat'));

[a1Len,e1Len,s1Len,f1Len,o1Len,c1Len,t1Len] = size(p.parameterCombinations);
[a2Len,e2Len,s2Len,f2Len,o2Len,c2Len,t2Len] = size(p.parameterCombinations2);

aLen = max(a1Len,a2Len);
eLen = max(e1Len,e2Len);
sLen = max(s1Len,s2Len);
fLen = max(f1Len,f2Len);
oLen = max(o1Len,o2Len);
cLen = max(c1Len,c2Len);
tLen = max(t1Len,t2Len);

% Get the unique values
aValsUnique = getVals(p.aValsUnique,p.aValsUnique2);
eValsUnique = getVals(p.eValsUnique,p.eValsUnique2);
sValsUnique = getVals(p.sValsUnique,p.sValsUnique2);
fValsUnique = getVals(p.fValsUnique,p.fValsUnique2);
oValsUnique = getVals(p.oValsUnique,p.oValsUnique2);
cValsUnique = getVals(p.cValsUnique,p.cValsUnique2);
tValsUnique = getVals(p.tValsUnique,p.tValsUnique2);


for a=1:aLen
    for e=1:eLen
        for s=1:sLen
            for f=1:fLen
                for o=1:oLen
                    for c=1:cLen
                        for t=1:tLen
                            
                            a1 = min(a,a1Len); a2 = min(a,a2Len);
                            e1 = min(e,e1Len); e2 = min(e,e2Len);
                            s1 = min(s,s1Len); s2 = min(s,s2Len);
                            f1 = min(f,f1Len); f2 = min(f,f2Len);
                            o1 = min(o,o1Len); o2 = min(o,o2Len);
                            c1 = min(c,c1Len); c2 = min(c,c2Len);
                            t1 = min(t,t1Len); t2 = min(t,t2Len);
                            
                            parameterCombinations{a,e,s,f,o,c,t} = intersect(p.parameterCombinations{a1,e1,s1,f1,o1,c1,t1},p.parameterCombinations2{a2,e2,s2,f2,o2,c2,t2});
                        end
                    end
                end
            end
        end
    end
end
end

function valsOut = getVals(valsUnique1,valsUnique2)

n1 = length(valsUnique1); n2 = length(valsUnique2);
if n2==1
    valsOut = valsUnique1;
elseif n1==1
    valsOut = valsUnique2;
else
    error('One of the two entries must be 1');
end
end