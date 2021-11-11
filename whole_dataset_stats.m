%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whole Dataset Stats:
% 
% Syntax:
% 
% 
% Inputs:
% set results path:
clear

dbPath = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\QT_results';
results_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\QT_results';

% 
% 
% Output:
% 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021  Sina Dabiri
% sdabiri@emory.edu
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
%% 


folderNames = ls(fullfile(dbPath, '*patient*'));
[numPatients, ~] = size(folderNames);

QT_table_colNames = {'subject','record', 'MedianQTlc_IQR_Fattahi',...
    'MedianQTlc_IQR_Li'};
colToRead = {'MedianQTlc_IQR_Fattahi_table',...
            'MedianQTlc_IQR_Li_table'};

% TODO: use datastructure to collect QTs
QT.subject={}; QT.record={};
QT.MedianQTlc_IQR_Fattahi=[]; QT.MedianQTlc_IQR_Li=[];
% rows = numPatients;
% QT = cell(rows);

for fn=1:80

    recordNames = ls(fullfile(dbPath, folderNames(fn,:),'s0*.csv'));
    [numRecords, ~] = size(recordNames);
    
   
    for rn=1:numRecords
        [~, baseFileName, extension] = fileparts(recordNames(rn, :));
        inPath = fullfile(dbPath, folderNames(fn,:), [baseFileName,extension]);
        
        opts = detectImportOptions(inPath);
        opts.SelectedVariableNames = colToRead;
        QT_read = readtable(inPath, opts);

        QT.subject(end+1,1) = {folderNames(fn,:)};
        QT.record(end+1,1) = {baseFileName(1:8)};
        QT.MedianQTlc_IQR_Fattahi(end+1,1) = table2array(QT_read(1,1));
        QT.MedianQTlc_IQR_Li(end+1,1) = table2array(QT_read(1,2));

        
    end
end
disp(QT)
%% General Stat measurements

% Discriptive stats
mean_Fattahi = mean(QT.MedianQTlc_IQR_Fattahi);
mean_Li = mean(QT.MedianQTlc_IQR_Li);

trimMean_Fattahi = trimmean(QT.MedianQTlc_IQR_Fattahi, 10);
trimMean_Li = trimmean(QT.MedianQTlc_IQR_Li, 10);

median_Fattahi = median(QT.MedianQTlc_IQR_Fattahi);
median_Li = median(QT.MedianQTlc_IQR_Li);

mode_Fattahi = mode(QT.MedianQTlc_IQR_Fattahi);
mode_Li = mode(QT.MedianQTlc_IQR_Li);

range_Fattahi = range(QT.MedianQTlc_IQR_Fattahi);
range_Li = range(QT.MedianQTlc_IQR_Li);

std_Fattahi = std(QT.MedianQTlc_IQR_Fattahi);
std_Li = std(QT.MedianQTlc_IQR_Li);

% printing table
QT_discriptive_stats = {'Mean', 'Trimmed Mean', 'Median', 'Mode', 'Range',...
    'Standard Deviation'};
QTlc_Gaussian = [mean_Fattahi trimMean_Fattahi median_Fattahi mode_Fattahi ...
    range_Fattahi std_Fattahi]';
QTlc_Wavelet = [median_Li trimMean_Li median_Li mode_Li range_Li std_Li]';

[discriptive_stats] = table(QTlc_Wavelet, QTlc_Gaussian, 'RowNames',QT_discriptive_stats)

% Calculate Box plot stats

% iqr
iqr_Gassian = iqr(QT.MedianQTlc_IQR_Fattahi);
iqr_Wavelet = iqr(QT.MedianQTlc_IQR_Li);

%iqr lower
iqr_Gassian_lower = prctile(QT.MedianQTlc_IQR_Fattahi, 25)-1.5*iqr_Gassian;
iqr_Wavelet_lower = prctile(QT.MedianQTlc_IQR_Li, 25)-1.5*iqr_Wavelet;

% iqr upper
iqr_Gassian_upper = prctile(QT.MedianQTlc_IQR_Fattahi, 75)+1.5*iqr_Gassian;
iqr_Wavelet_upper = prctile(QT.MedianQTlc_IQR_Li, 75)+1.5*iqr_Wavelet;

% box plot stats table: relative statistics
QT_relative_statistics = {'max','3rd IQR', 'Median','1st IQR', 'min'};
QTlc_Gaussian_iqr = [max(QT.MedianQTlc_IQR_Fattahi) iqr_Gassian_upper...
    median_Fattahi iqr_Gassian_lower min(QT.MedianQTlc_IQR_Fattahi)]';
QTlc_Wavelet_iqr = [max(QT.MedianQTlc_IQR_Li) iqr_Wavelet_upper...
    median_Li iqr_Wavelet_lower min(QT.MedianQTlc_IQR_Li)]';
boxPlot_stats = table(QTlc_Wavelet_iqr, QTlc_Gaussian_iqr, 'RowNames',...
    QT_relative_statistics)

%% Box plot

boxplot([QT.MedianQTlc_IQR_Li, QT.MedianQTlc_IQR_Fattahi],...
    'Labels',{'Median QTlc Wavelet', 'Median QTlc Gaussian'});
ylabel('time (seconds)');
title('PTB database Median QTlc')
% outPath = 'stats.csv';