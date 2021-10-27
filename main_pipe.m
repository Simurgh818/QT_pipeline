%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QT Measurement Pipeline ECG signal
% 
% Syntax:
% 
% 
% Inputs:
% 
% 
% Output:
% 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021  Sina Dabiri
% sdabiri@emory.edu
% 
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

dbPath = 'C:/Users/sinad/wfdb/database/ptbdb/';
folderNames = ls(fullfile(dbPath, '*patient*'));
[numPatients, ~] = size(folderNames);
for fn=1:2
%     TODO: add inner loop, to loop through each subject's records
    recordNames = ls(fullfile(dbPath, folderNames(fn,:),'s0*.dat'));
    [numRecords, ~] = size(recordNames);
    for rn=1:numRecords
        [~, baseFileName, extension] = fileparts(recordNames(rn, :));
        inPath = ['ptbdb/', folderNames(fn,:),'/', baseFileName];

        outPath = 'preprocessed.csv';
        preprocessing(inPath, outPath);

        processedPath = outPath; 
        [fPath,fName,fExt]=fileparts(inPath);
        figPath = fullfile('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\QT_results', folderNames(fn,:));
        if ~exist(figPath, 'dir')
            mkdir(figPath)
        end

        QT_measurements(processedPath, fName, figPath);
        
    end
end



