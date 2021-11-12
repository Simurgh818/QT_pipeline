function [QT] = whole_dataset_stats(results_path)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whole Dataset Stats:
% 
% Syntax: whole_dataset_stats
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
% ToDo: save the results in a CSV in the parent result folder.

folderNames = ls(fullfile(dbPath, '*patient*'));
[numPatients, ~] = size(folderNames);


colToRead = {'MedianQTlc_IQR_Fattahi_table',...
            'MedianQTlc_IQR_Li_table'};

% TODO: use datastructure to collect QTs
QT.subject={}; QT.record={};
QT.MedianQTlc_IQR_Fattahi=[]; QT.MedianQTlc_IQR_Li=[];
% rows = numPatients;
% QT = cell(rows);

for fn=1:numPatients

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
%% Stat measurements

% Chi-squared Test

% chi-squared probability density function distribution
[Gaussian_size ,~] = size(QT.MedianQTlc_IQR_Fattahi);
x = 1:Gaussian_size;
y = chi2pdf(x,1);

% Normal distribution check: 
[h, p] = chi2gof(QT.MedianQTlc_IQR_Fattahi,'Alpha',0.01);
[h2, p2] = chi2gof(QT.MedianQTlc_IQR_Li, 'Alpha',0.01);
if h==1 && h2==1
    disp('reject null hypothesis, not a normal distribution');
    fprintf('The Gaussian p-value is %d, and the Wavelet p-value is %d.',...
        p, p2);
else
    disp('null hypotheiss not reject, a normal distribution.');
end

% chi-squared test:
% chi_sqrd = sum((QT.MedianQTlc_IQR_Fattahi-QT.MedianQTlc_IQR_Li).^2/...
%     QT.MedianQTlc_IQR_Fattahi,[1, Gaussian_size]);
% if chi_sqrd == 0
%     disp('No difference between Gaussian and Wavelet')
% else
%     disp('There is a difference by chi-squared method')
% end

%% plot

% figure(1)
% scatter(x, QT.MedianQTlc_IQR_Fattahi, 'r', 'filled');
% ylabel('Second')
% title('The Gaussian Model distribution for Median QTlc IQR')
% 
% figure(2)
% x2 = size(QT.MedianQTlc_IQR_Li);
% scatter(1:x2, QT.MedianQTlc_IQR_Li, 'b', 'filled');
% ylabel('Second')
% title('The Wavelet method distribution for Median QTlc IQR')

figure(3)
plot(x,y);
xlabel('observations')
ylabel('probability density')
%% Save results for the whole dataset in a CSV


QT_table_colNames = {'subject','record', 'MedianQTlc_IQR_Fattahi',...
    'MedianQTlc_IQR_Li'};
QT_table = table(QT.subject, QT.record, QT.MedianQTlc_IQR_Fattahi, ...
    QT.MedianQTlc_IQR_Li, 'VariableNames',QT_table_colNames);


csv_fileName = 'QT_stats.csv';
fileName = fullfile(results_path, csv_fileName);
writetable(QT_table, fileName);
end