%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QT Measurement PipeLie ECG signal
% 
% Syntax: main_pipe
% 
% 
% Inputs:
% 
% set database path:
clear
close all;

dbPath = 'C:\Users\sinad\wfdb\10.6.2\database\qtdb\physionet.org\files\qtdb';
nChannels = 2; % number of channels to read
% 
% set results path:
results_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\QT_results\QTdb_results';
 
% Output: It saves scatter plot for each record in a subfolder of each 
%         subject in the results_path. It also saves a CSV with the
%         following measurements for Gaussian model correct QT: median, mean, median IQR,
%         mean RR and median RR. Also, for the wavelet method correct QT:
%         median, median IQR, gauss jQRS, Median RR, jQRS RR.
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
%
%       Method_1: Fattahi Gaussian model
%       Method_2: Li wavelet method  
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

% ToDo: Need to generalize a way to list patient/subject subfolders for the
% case there is subfolder and the case there is no subfolders: use try and
% catch



% folderNames = ls(fullfile(dbPath));
dirNames = dir(fullfile(dbPath));
folderNames = dirNames([dirNames(:).isdir]);
folderNames = folderNames(~ismember({folderNames(:).name},{'.','..'}));

[numPatients, ~] = size(folderNames);


for fn=1:numPatients
    fprintf('Currently Processing subject: %s \n', folderNames(fn,:).name);
    recordNames = dir(fullfile(dbPath, folderNames(fn,:).name,'*.dat'));

    [numRecords, ~] = size(recordNames);
    for rn=1:1
        fprintf('Currently Processing record #: %s \n', recordNames(rn,:).name);
        [~, baseFileName, extension] = fileparts(recordNames(rn, :).name);
        pathSplit = split(dbPath, '\');
        inPath = [pathSplit{end, 1}, '/', folderNames(fn,:).name,'/', baseFileName];

        outPath = 'preprocessed.csv';
        [fs] = preprocessing(inPath, outPath, nChannels);

        processedPath = outPath; 
        [fPath,fName,fExt]=fileparts(inPath);
        figPath = fullfile(results_path, folderNames(fn,:).name);
        if ~exist(figPath, 'dir')
            mkdir(figPath)
        end
%         Reading human annotations
        [humanQT] = humanAnnotations(inPath,'q1c', fName, figPath, fs, nChannels);
%       Method_1: Fattahi
%       Method_2: Li    
        [Method_2, Method_1] = QT_measurements(outPath, fName, figPath, nChannels, fs, humanQT);



    end
end

% catch
%     if isempty(folderNames)
%         fprintf('Dataset did not load properly. Please check the path.');
%     end
% 
% end

% QT = whole_dataset_stats(results_path,nChannels);


