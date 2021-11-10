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

for fn=1:100

    recordNames = ls(fullfile(dbPath, folderNames(fn,:),'s0*.csv'));
    [numRecords, ~] = size(recordNames);
    
   
    for rn=1:numRecords
        [~, baseFileName, extension] = fileparts(recordNames(rn, :));
        inPath = fullfile(dbPath, folderNames(fn,:), [baseFileName,extension]);
        
        opts = detectImportOptions(inPath);
        opts.SelectedVariableNames = colToRead;
        QT_read = readtable(inPath, opts);
%         TODO: use a structure instead

%         QT_table = table({folderNames(fn,:)}, {baseFileName},...
%             table2cell(QT(1,1)), table2cell(QT(1,2)),...
%             'VariableNames', QT_table_colNames);
        row = fn*rn;
        QT.subject(end+1,1) = {folderNames(fn,:)};
        QT.record(end+1,1) = {baseFileName(1:8)};
        QT.MedianQTlc_IQR_Fattahi(end+1,1) = table2array(QT_read(1,1));
        QT.MedianQTlc_IQR_Li(end+1,1) = table2array(QT_read(1,2));

%         QT_m = table2array(QT(1,:));
%         QT_record = [folderNames; rn; ]
        
    end
end
disp(QT)
%% Box plot
boxplot([QT.MedianQTlc_IQR_Fattahi, QT.MedianQTlc_IQR_Li])
xlabel(['Median QTlc Fattahi ', ' Median QTlc Li']);
ylabel('time (seconds)');
title('PTB database Median QTlc')
% outPath = 'stats.csv';