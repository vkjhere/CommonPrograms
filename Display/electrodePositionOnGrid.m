function [row,column,electrodeArray] = electrodePositionOnGrid(electrodeNum,gridType,subjectName)

if ~exist('subjectName','var');         subjectName=[];                 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(gridType,'EEG')
    electrodeArray = ...
         [00 00 00 00 01 61 02 00 00 00 00;
         00 00 53 00 39 00 40 00 54 00 00;
         00 11 47 03 33 17 34 04 48 12 00;
         29 55 25 41 21 20 22 42 26 56 30;
         00 13 49 05 35 18 36 06 50 14 00;
         31 57 27 43 23 62 24 44 28 58 32;
         00 15 51 07 37 19 38 08 52 16 00;
         00 00 59 00 45 63 46 00 60 00 00;
         00 00 00 00 09 64 10 00 00 00 00;
         00 00 00 00 00 00 00 00 00 00 00];
     
     
     if electrodeNum<1 || electrodeNum>96
         disp('Electrode Number out of range');
         row=1;column=1;
     else
         [row,column] = find(electrodeArray==electrodeNum);
     end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Microelectrode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if strcmpi(gridType,'Microelectrode')
    if strcmp(subjectName,'abu') || strcmp(subjectName,'rafiki') || isempty(subjectName)
        
        electrodeArray = ...
            [00 02 01 03 04 06 08 10 14 00;
            65 66 33 34 07 09 11 12 16 18;
            67 68 35 36 05 17 13 23 20 22;
            69 70 37 38 48 15 19 25 27 24;
            71 72 39 40 42 50 54 21 29 26;
            73 74 41 43 44 46 52 62 31 28;
            75 76 45 47 51 56 58 60 64 30;
            77 78 82 49 53 55 57 59 61 32;
            79 80 84 86 87 89 91 94 63 95;
            00 81 83 85 88 90 92 93 96 00];
        
    elseif strcmp(subjectName,'alpa')
        
        electrodeArray = ...
           [00 88 78 68 58 48 38 28 18 00;
            96 87 77 67 57 47 37 27 17 08;
            95 86 76 66 56 46 36 26 16 07;
            94 85 75 65 55 45 35 25 15 06;
            93 84 74 64 54 44 34 24 14 05;
            92 83 73 63 53 43 33 23 13 04;
            91 82 72 62 52 42 32 22 12 03;
            90 81 71 61 51 41 31 21 11 02;
            89 80 70 60 50 40 30 20 10 01;
            00 79 69 59 49 39 29 19 09 00];
    end
    
    if electrodeNum<1 || electrodeNum>96
        disp('Electrode Number out of range');
        row=1;column=1;
    else
        [row,column] = find(electrodeArray==electrodeNum);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ECoG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(gridType,'ECoG')
    
    goodElectrodeList = [1:34 65:96];
    electrodeArray = ...
        [00 89 81 73 65 25 17 09 01 00;
        00 90 82 74 66 26 18 10 02 00;
        00 91 83 75 67 27 19 11 03 00;
        34 92 84 76 68 28 20 12 04 33;
        00 93 85 77 69 29 21 13 05 00;
        00 94 86 78 70 30 22 14 06 00;
        00 95 87 79 71 31 23 15 07 00;
        00 96 88 80 72 32 24 16 08 00];
    
    if isempty(setdiff(goodElectrodeList,electrodeNum))
        disp('Electrode Number out of range');
    else
        [row,column] = find(electrodeArray==electrodeNum);
    end
end
end