% This program generates the parameterCombinations variable from the
% stimResults
function parameterCombinations = getDisplayCombinationsGRF(folderOut,goodStimNums,radiusToSigmaRatio)

if ~exist('radiusToSigmaRatio','var');      radiusToSigmaRatio=3;       end

load(fullfile(folderOut,'stimResults.mat'));

% Five parameters are chosen:
% 1. Azimuth
% 2. Elevation
% 3. Sigma, Radius 
% 4. Spatial Frequency
% 5. Orientation
% 6. Contrast
% 7. Temporal Frequency

% Parameters index
parameters{1} = 'azimuth';
parameters{2} = 'elevation';
parameters{3} = 'sigma';
parameters{4} = 'spatialFrequency';
parameters{5} = 'orientation';
parameters{6} = 'contrast';
parameters{7} = 'temporalFrequency'; %#ok<NASGU>

if stimResults.side == 3 % Plaids
    aValsAll  = stimResults.azimuth(1:2:end);
    eValsAll  = stimResults.elevation(1:2:end);
    if isfield(stimResults,'radius')
        sValsAll  = stimResults.radius(1:2:end)/radiusToSigmaRatio;
    else
        sValsAll  = stimResults.sigma(1:2:end);
    end

    fValsAll  = stimResults.spatialFrequency(1:2:end);
    oValsAll  = stimResults.orientation(1:2:end);
    cValsAll  = stimResults.contrast(1:2:end);
    tValsAll  = stimResults.temporalFrequency(1:2:end);

    aValsAll2 = stimResults.azimuth(2:2:end);
    eValsAll2 = stimResults.elevation(2:2:end);
    if isfield(stimResults,'radius')
        sValsAll2  = stimResults.radius(2:2:end)/radiusToSigmaRatio;
    else
        sValsAll2  = stimResults.sigma(2:2:end);
    end

    fValsAll2 = stimResults.spatialFrequency(2:2:end);
    oValsAll2 = stimResults.orientation(2:2:end);
    cValsAll2 = stimResults.contrast(2:2:end);
    tValsAll2 = stimResults.temporalFrequency(2:2:end);

else
    aValsAll  = stimResults.azimuth;
    eValsAll  = stimResults.elevation;
    if isfield(stimResults,'radius')
        sValsAll  = stimResults.radius/radiusToSigmaRatio;
    else
        sValsAll  = stimResults.sigma;
    end

    fValsAll  = stimResults.spatialFrequency;
    oValsAll  = stimResults.orientation;
    cValsAll  = stimResults.contrast;
    tValsAll  = stimResults.temporalFrequency;
end

if ~isempty(aValsAll)

    aValsGood = aValsAll(goodStimNums);
    eValsGood = eValsAll(goodStimNums);
    sValsGood = sValsAll(goodStimNums);
    fValsGood = fValsAll(goodStimNums);
    oValsGood = oValsAll(goodStimNums);
    cValsGood = cValsAll(goodStimNums);
    tValsGood = tValsAll(goodStimNums);

    aValsUnique = unique(aValsGood); aLen = length(aValsUnique);
    eValsUnique = unique(eValsGood); eLen = length(eValsUnique);
    sValsUnique = unique(sValsGood); sLen = length(sValsUnique);
    fValsUnique = unique(fValsGood); fLen = length(fValsUnique);
    oValsUnique = unique(oValsGood); oLen = length(oValsUnique);
    cValsUnique = unique(cValsGood); cLen = length(cValsUnique);
    tValsUnique = unique(tValsGood); tLen = length(tValsUnique);

    % display
    if stimResults.side == 3 % Plaids
        disp('Left plaid component');
    end
    disp(['Number of unique azimuths: ' num2str(aLen)]);
    disp(['Number of unique elevations: ' num2str(eLen)]);
    disp(['Number of unique sigmas: ' num2str(sLen)]);
    disp(['Number of unique Spatial freqs: ' num2str(fLen)]);
    disp(['Number of unique orientations: ' num2str(oLen)]);
    disp(['Number of unique contrasts: ' num2str(cLen)]);
    disp(['Number of unique temporal freqs: ' num2str(tLen)]);

    % If more than one value, make another entry with all values
    if (aLen > 1);           aLen=aLen+1;                    end
    if (eLen > 1);           eLen=eLen+1;                    end
    if (sLen > 1);           sLen=sLen+1;                    end
    if (fLen > 1);           fLen=fLen+1;                    end
    if (oLen > 1);           oLen=oLen+1;                    end
    if (cLen > 1);           cLen=cLen+1;                    end
    if (tLen > 1);           tLen=tLen+1;                    end

    allPos = 1:length(goodStimNums);
    disp(['total combinations: ' num2str((aLen)*(eLen)*(sLen)*(fLen)*(oLen)*(cLen)*(tLen))]);

    for a=1:aLen
        if a==aLen
            aPos = allPos;
        else
            aPos = find(aValsGood == aValsUnique(a));
        end

        for e=1:eLen
            if e==eLen
                ePos = allPos;
            else
                ePos = find(eValsGood == eValsUnique(e));
            end

            for s=1:sLen
                if s==sLen
                    sPos = allPos;
                else
                    sPos = find(sValsGood == sValsUnique(s));
                end

                for f=1:fLen
                    if f==fLen
                        fPos = allPos;
                    else
                        fPos = find(fValsGood == fValsUnique(f));
                    end

                    for o=1:oLen
                        if o==oLen
                            oPos = allPos;
                        else
                            oPos = find(oValsGood == oValsUnique(o));
                        end
                        
                        for c=1:cLen
                            if c==cLen
                                cPos = allPos;
                            else
                                cPos = find(cValsGood == cValsUnique(c));
                            end
                            
                            for t=1:tLen
                                if t==tLen
                                    tPos = allPos;
                                else
                                    tPos = find(tValsGood == tValsUnique(t));
                                end


                                aePos = intersect(aPos,ePos);
                                aesPos = intersect(aePos,sPos);
                                aesfPos = intersect(aesPos,fPos);
                                aesfoPos = intersect(aesfPos,oPos);
                                aesfocPos = intersect(aesfoPos,cPos);
                                aesfoctPos = intersect(aesfocPos,tPos);
                                parameterCombinations{a,e,s,f,o,c,t} = aesfoctPos; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end
    end

    % save
    save(fullfile(folderOut,'parameterCombinations.mat'),'parameters','parameterCombinations', ...
        'aValsUnique','eValsUnique','sValsUnique','fValsUnique','oValsUnique','cValsUnique','tValsUnique');

    if stimResults.side == 3 % Plaids
        
        aValsGood = aValsAll2(goodStimNums);
        eValsGood = eValsAll2(goodStimNums);
        sValsGood = sValsAll2(goodStimNums);
        fValsGood = fValsAll2(goodStimNums);
        oValsGood = oValsAll2(goodStimNums);
        cValsGood = cValsAll2(goodStimNums);
        tValsGood = tValsAll2(goodStimNums);

        aValsUnique2 = unique(aValsGood); aLen = length(aValsUnique2);
        eValsUnique2 = unique(eValsGood); eLen = length(eValsUnique2);
        sValsUnique2 = unique(sValsGood); sLen = length(sValsUnique2);
        fValsUnique2 = unique(fValsGood); fLen = length(fValsUnique2);
        oValsUnique2 = unique(oValsGood); oLen = length(oValsUnique2);
        cValsUnique2 = unique(cValsGood); cLen = length(cValsUnique2);
        tValsUnique2 = unique(tValsGood); tLen = length(tValsUnique2);

        % display
        disp('Right plaid component');
        disp(['Number of unique azimuths: ' num2str(aLen)]);
        disp(['Number of unique elevations: ' num2str(eLen)]);
        disp(['Number of unique sigmas: ' num2str(sLen)]);
        disp(['Number of unique Spatial freqs: ' num2str(fLen)]);
        disp(['Number of unique orientations: ' num2str(oLen)]);
        disp(['Number of unique contrasts: ' num2str(cLen)]);
        disp(['Number of unique temporal freqs: ' num2str(tLen)]);

        % If more than one value, make another entry with all values
        if (aLen > 1);           aLen=aLen+1;                    end
        if (eLen > 1);           eLen=eLen+1;                    end
        if (sLen > 1);           sLen=sLen+1;                    end
        if (fLen > 1);           fLen=fLen+1;                    end
        if (oLen > 1);           oLen=oLen+1;                    end
        if (cLen > 1);           cLen=cLen+1;                    end
        if (tLen > 1);           tLen=tLen+1;                    end

        allPos = 1:length(goodStimNums);
        disp(['total combinations: ' num2str((aLen)*(eLen)*(sLen)*(fLen)*(oLen)*(cLen)*(tLen))]);

        for a=1:aLen
            if a==aLen
                aPos = allPos;
            else
                aPos = find(aValsGood == aValsUnique2(a));
            end

            for e=1:eLen
                if e==eLen
                    ePos = allPos;
                else
                    ePos = find(eValsGood == eValsUnique2(e));
                end

                for s=1:sLen
                    if s==sLen
                        sPos = allPos;
                    else
                        sPos = find(sValsGood == sValsUnique2(s));
                    end

                    for f=1:fLen
                        if f==fLen
                            fPos = allPos;
                        else
                            fPos = find(fValsGood == fValsUnique2(f));
                        end

                        for o=1:oLen
                            if o==oLen
                                oPos = allPos;
                            else
                                oPos = find(oValsGood == oValsUnique2(o));
                            end

                            for c=1:cLen
                                if c==cLen
                                    cPos = allPos;
                                else
                                    cPos = find(cValsGood == cValsUnique2(c));
                                end

                                for t=1:tLen
                                    if t==tLen
                                        tPos = allPos;
                                    else
                                        tPos = find(tValsGood == tValsUnique2(t));
                                    end


                                    aePos = intersect(aPos,ePos);
                                    aesPos = intersect(aePos,sPos);
                                    aesfPos = intersect(aesPos,fPos);
                                    aesfoPos = intersect(aesfPos,oPos);
                                    aesfocPos = intersect(aesfoPos,cPos);
                                    aesfoctPos = intersect(aesfocPos,tPos);
                                    parameterCombinations2{a,e,s,f,o,c,t} = aesfoctPos; %#ok<NASGU,AGROW>
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % save
        save(fullfile(folderOut,'parameterCombinations.mat'),'parameterCombinations2', ...
            'aValsUnique2','eValsUnique2','sValsUnique2','fValsUnique2','oValsUnique2','cValsUnique2','tValsUnique2', ...
            '-append');
    
    end
    
end
end