%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QT Measurement PipeLie ECG signal
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
% Add OSET's Tools folder to the Matlab path
% addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\OSET\Tools\')
%
% 1- Non-model method: Dr. Li wavelet method
%  https://github.com/cliffordlab/QTestimation.git
% Add Dr. Li's QT estimator folder to Matlab path
% addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\QTestimation\QTestimation\QT_for_Alivecor\')
%
% 2- Add Mr. Fattahi's QT folder to Matlab path
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval
% addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\OSET\UnderDevelopment\QTinterval');
%
% set database path:
dbPath = 'C:/Users/sinad/wfdb/10.6.2/database/ptbdb/';
nChannels = 14; % number of channels to read
%
% set results path:
results_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\QT_results';
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

folderNames = ls(fullfile(dbPath, '*patient*'));
[numPatients, ~] = size(folderNames);
for fn=105:numPatients
    fprintf('Currently Processing subject: %s \n', folderNames(fn,:));
    recordNames = ls(fullfile(dbPath, folderNames(fn,:),'s0*.dat'));
    [numRecords, ~] = size(recordNames);
    for rn=1:numRecords
        [~, baseFileName, extension] = fileparts(recordNames(rn, :));
        inPath = ['ptbdb/', folderNames(fn,:),'/', baseFileName];

        outPath = 'preprocessed.csv';
        [fs] = preprocessing(inPath, outPath, nChannels);

        processedPath = outPath; 
        [fPath,fName,fExt]=fileparts(inPath);
        figPath = fullfile(results_path, folderNames(fn,:));
        if ~exist(figPath, 'dir')
            mkdir(figPath)
        end

        [Li, Fattahi] = QT_measurements(processedPath, fName, figPath, nChannels, fs);
        
    end
end

% TODO: add a whole dataset boxplot

