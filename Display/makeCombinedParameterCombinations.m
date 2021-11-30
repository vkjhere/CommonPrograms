% When both Gabor0 and Gabor1 stimuli are used together, say in case of a
% plaid, it is possible to combine the parameterCombinations when only one
% out of two gabor paramaters takes more than 1 value. For example, if grating
% 1 has 4 orientations and 1 TF while grating2 has 1 orientation and 10
% TFs, this program will return a combined parameterCombinations matrix
% with 4 orientations and 10 TFs.

% Adding optional parameter that specifies which of the two sides to use.
% This should be of the form sideChoice.a/e/s/f/o/c/t, each taking the
% value 1 or 2.

function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = makeCombinedParameterCombinations(folderExtract,sideChoice)

if ~exist('sideChoice','var');      sideChoice=[];                      end

p = load(fullfile(folderExtract,'parameterCombinations.mat'));
[a1Len,e1Len,s1Len,f1Len,o1Len,c1Len,t1Len] = size(p.parameterCombinations); % Side 1
[a2Len,e2Len,s2Len,f2Len,o2Len,c2Len,t2Len] = size(p.parameterCombinations2); % Side 2

if isempty(sideChoice)
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
    
else
    if sideChoice.a==1
        aValsUnique = p.aValsUnique; aLen = a1Len;
    elseif sideChoice.a==2
        aValsUnique = p.aValsUnique2; aLen = a2Len;
    else
        error('Choice parameter should be 1 or 2');
    end
    
    if sideChoice.e==1
        eValsUnique = p.eValsUnique; eLen = e1Len;
    elseif sideChoice.e==2
        eValsUnique = p.eValsUnique2; eLen = e2Len;
    else
        error('Choice parameter should be 1 or 2');
    end
    
    if sideChoice.s==1
        sValsUnique = p.sValsUnique; sLen = s1Len;
    elseif sideChoice.s==2
        sValsUnique = p.sValsUnique2; sLen = s2Len;
    else
        error('Choice parameter should be 1 or 2');
    end
    
    if sideChoice.f==1
        fValsUnique = p.fValsUnique; fLen = f1Len;
    elseif sideChoice.f==2
        fValsUnique = p.fValsUnique2; fLen = f2Len;
    else
        error('Choice parameter should be 1 or 2');
    end
    
    if sideChoice.o==1
        oValsUnique = p.oValsUnique; oLen = o1Len;
    elseif sideChoice.o==2
        oValsUnique = p.oValsUnique2; oLen = o2Len;
    else
        error('Choice parameter should be 1 or 2');
    end
    
    if sideChoice.c==1
        cValsUnique = p.cValsUnique; cLen = c1Len;
    elseif sideChoice.c==2
        cValsUnique = p.cValsUnique2; cLen = c2Len;
    else
        error('Choice parameter should be 1 or 2');
    end
    
    if sideChoice.t==1
        tValsUnique = p.tValsUnique; tLen = t1Len;
    elseif sideChoice.t==2
        tValsUnique = p.tValsUnique2; tLen = t2Len;
    else
        error('Choice parameter should be 1 or 2');
    end
    
    for a=1:aLen
        for e=1:eLen
            for s=1:sLen
                for f=1:fLen
                    for o=1:oLen
                        for c=1:cLen
                            for t=1:tLen
                                
                                if sideChoice.a==1
                                    a1 = a; a2 = a2Len;
                                else
                                    a1 = a1Len; a2 = a;
                                end
                                
                                if sideChoice.e==1
                                    e1 = e; e2 = e2Len;
                                else
                                    e1 = e1Len; e2 = e;
                                end
                                
                                if sideChoice.s==1
                                    s1 = s; s2 = s2Len;
                                else
                                    s1 = s1Len; s2 = s;
                                end
                                
                                if sideChoice.f==1
                                    f1 = f; f2 = f2Len;
                                else
                                    f1 = f1Len; f2 = f;
                                end
                                
                                if sideChoice.o==1
                                    o1 = o; o2 = o2Len;
                                else
                                    o1 = o1Len; o2 = o;
                                end
                                
                                if sideChoice.c==1
                                    c1 = c; c2 = c2Len;
                                else
                                    c1 = c1Len; c2 = c;
                                end
                                
                                if sideChoice.t==1
                                    t1 = t; t2 = t2Len;
                                else
                                    t1 = t1Len; t2 = t;
                                end
                                
                                parameterCombinations{a,e,s,f,o,c,t} = intersect(p.parameterCombinations{a1,e1,s1,f1,o1,c1,t1},p.parameterCombinations2{a2,e2,s2,f2,o2,c2,t2}); %#ok<*AGROW>
                            end
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