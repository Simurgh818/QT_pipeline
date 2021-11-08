%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whole Dataset Stats:
% 
% Syntax:
% 
% 
% Inputs:
% set results path:
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
    'MedianQTlc_IQR_Lin'};
colToRead = {'MedianQTlc_IQR_Fattahi_table',...
            'MedianQTlc_IQR_Lin_table'};

for fn=1:numPatients

    recordNames = ls(fullfile(dbPath, folderNames(fn,:),'s0*.csv'));
    [numRecords, ~] = size(recordNames);
    for rn=1:numRecords
        [~, baseFileName, extension] = fileparts(recordNames(rn, :));
        inPath = fullfile(dbPath, folderNames(fn,:), [baseFileName,extension]);
        
        opts = detectImportOptions(inPath);
        opts.SelectedVariableNames = colToRead;
        QT = readtable(inPath, opts);
%         TODO: use a structure instead

        QT_table = table({folderNames(fn,:)}, {baseFileName},...
            table2cell(QT(1,1)), table2cell(QT(1,2)),...
            'VariableNames', QT_table_colNames);
%         QT_m = table2array(QT(1,:));
%         QT_record = [folderNames; rn; ]
        
    end
end
disp(QT_table)

% outPath = 'stats.csv';