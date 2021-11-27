function [QT] = whole_dataset_stats(results_path, numChannels)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whole Dataset Stats:
% 
% Syntax: whole_dataset_stats
% 
% 
% Inputs:
% set results path:
% dbPath = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\QT_results\QTdb_results';
% results_path = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\QT_results';

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

dirNames = dir(fullfile(results_path));
folderNames = dirNames([dirNames(:).isdir]);
folderNames = folderNames(~ismember({folderNames(:).name},{'.','..'}));
[numPatients, ~] = size(folderNames);


colToRead = {'RR_mean','RR_median', 'QTc1_median', 'QTc1_mean',...
    'QTc1_median_IQR_Method_1_table', 'RR_median_wavelet', 'RR_mean_wavelet', 'RR_jQRS',...
    'QTc1_median_wavelet', 'QTc1_median_wavelet_SQI', 'QTc1_mean_jQRS',...
    'GaussQTlc_jQRS', 'GaussQTlc_SQI', 'QTc1_median_IQR_Method_2_table'};

% TODO: use datastructure to collect QTs
QT.subject={}; QT.record={}; QT.channelNum=[] ; QT.RR_mean_method1=[];QT.RR_median_method1=[];
QT.QTc1_median_method1=[]; QT.QTc1_mean_method1=[];
QT.QTc1_median_IQR_Method_1=[]; 
QT.RR_median_method2=[]; QT.RR_jQRS_method2=[]; QT.RR_mean_method2=[];
QT.QTc1_median_method2=[]; QT.QTc1_median_method2_SQI=[]; QT.QTc1_mean_method2_jQRS=[];
QT.QTc1_jQRS_method2_Gauss=[]; QT.QTc1_jQRS_method2_Gauss_SQI=[];
QT.QTc1_median_IQR_Method_2=[]; QT.RR_median_human=[]; QT.RR_mean_human=[];
QT.QTc1_median_human=[]; QT.QTc1_mean_human=[]; QT.RMSE_method_1=[]; 
QT.RMSE_method_2=[]; 
% rows = numPatients;
% QT = cell(rows);

% Human QT csv
colToReadHuman = {'RR_median_human', 'RR_mean_human',...
    'QTc1_median_human', 'QTc1_mean_human'};

for fn=1:numPatients
    recordNames = dir(fullfile(results_path, folderNames(fn,:).name,'*QT_results.csv'));
    recordNamesHuman = dir(fullfile(results_path, folderNames(fn,:).name,'*humanQT.csv'));
    [numRecords, ~] = size(recordNames);
    
   
    for rn=1:numRecords
        [~, baseFileName, extension] = fileparts(recordNames(rn, :).name);
        inPath = fullfile(results_path, folderNames(fn,:).name, [baseFileName,extension]);
        [~, baseFileNameHuman, ~] = fileparts(recordNamesHuman(rn, :).name);
        inPathHuman = fullfile(results_path,folderNames(fn,:).name, [baseFileNameHuman,extension]);

        opts = detectImportOptions(inPath);
        opts.SelectedVariableNames = colToRead;
        QT_read = readtable(inPath, opts);
        
        optsHuman = detectImportOptions(inPathHuman);
        optsHuman.SelectedVariableNames = colToReadHuman;
        QT_read_human = readtable(inPathHuman, optsHuman);

        for ch=1:numChannels
    
            QT.subject(end+1,1) = {folderNames(fn,:).name};
            split_bfn = split(baseFileName,'_');
            QT.record(end+1,1) = split_bfn(1,:);
            QT.channelNum(end+1,1) = ch;
            QT.RR_mean_method1(end+1,1) = table2array(QT_read(ch,1));
            QT.RR_median_method1(end+1,1) = table2array(QT_read(ch,2));
            QT.QTc1_median_method1(end+1,1) = table2array(QT_read(ch,3));
            QT.QTc1_mean_method1(end+1,1) = table2array(QT_read(ch,4));
            QT.QTc1_median_IQR_Method_1(end+1,1) = table2array(QT_read(ch,5));
            QT.RR_median_method2(end+1,1) = table2array(QT_read(ch,6));
            QT.RR_mean_method2(end+1,1) = table2array(QT_read(ch,7));
            QT.RR_jQRS_method2(end+1,1)= table2array(QT_read(ch,8));

            QT.QTc1_median_method2(end+1,1) = table2array(QT_read(ch,9));
            QT.QTc1_median_method2_SQI(end+1,1) = table2array(QT_read(ch,10));
            QT.QTc1_mean_method2_jQRS(end+1,1) = table2array (QT_read(ch, 11));
            QT.QTc1_jQRS_method2_Gauss(end+1,1) = table2array(QT_read(ch, 12));
            QT.QTc1_jQRS_method2_Gauss_SQI(end+1,1) = table2array(QT_read(ch,13));
            QT.QTc1_median_IQR_Method_2(end+1,1) = table2array(QT_read(ch,14));
            
            QT.RR_median_human(end+1,1) = table2array(QT_read_human(ch,1));
            QT.RR_mean_human(end+1,1) = table2array(QT_read_human(ch, 2));
            QT.QTc1_median_human(end+1,1) = table2array(QT_read_human(ch, 3));
            QT.QTc1_mean_human(end+1,1) = table2array(QT_read_human(ch, 4));

        end
    end
end
%% Stat measurements

% Normalized Root Mean Squared Error for Method 1 and 2 vs. Human
% annotation
% fit_nrmse = goodnessOfFit(yc,yrefc,'NRMSE')
% QT.RMSE_method_1 = sqrt(nanmean((QT.QTc1_median_IQR_Method_1-QT.QTc1_median_human).^2));
% QT.RMSE_method_2 = sqrt(nanmean((QT.QTc1_median_IQR_Method_2-QT.QTc1_median_human).^2));
% fprintf('The RMSE between human QT and Method 1 is %1.3g, and for Method 2 is %1.3g.\n',...
%         QT.RMSE_method_1, QT.RMSE_method_2);

% channel 1
QTc1_median_method_1_ch1=QT.QTc1_median_method1(QT.channelNum==1);
QTc1_median_method_2_ch1=QT.QTc1_median_method2(QT.channelNum==1);
QTc1_median_human_ch1= QT.QTc1_median_human(QT.channelNum==1);
% channel 2
QTc1_median_method_1_ch2=QT.QTc1_median_method1(QT.channelNum==2);
QTc1_median_method_2_ch2=QT.QTc1_median_method2(QT.channelNum==2);
QTc1_median_human_ch2= QT.QTc1_median_human(QT.channelNum==2);

QT.RMSE_method_1(1,1) = sqrt(nanmean((QTc1_median_method_1_ch1-QTc1_median_human_ch1).^2));
QT.RMSE_method_2(1,1) = sqrt(nanmean((QTc1_median_method_2_ch1-QTc1_median_human_ch1).^2));
QT.RMSE_method_1(2,1) = sqrt(nanmean((QTc1_median_method_1_ch2-QTc1_median_human_ch2).^2));
QT.RMSE_method_2(2,1) = sqrt(nanmean((QTc1_median_method_2_ch2-QTc1_median_human_ch2).^2));
disp(QT)
fprintf('The RMSE between lead 1 human QTc1 and Method 1 is %1.3g, and for Method 2 is %1.3g.\n',...
        QT.RMSE_method_1(1,1), QT.RMSE_method_2(1,1));
fprintf('The RMSE between lead 2 human QTc1 and Method 1 is %1.3g, and for Method 2 is %1.3g.\n',...
        QT.RMSE_method_1(2,1), QT.RMSE_method_2(2,1));

% Chi-squared Test

% chi-squared probability density function distribution
[Gaussian_size ,~] = size(QT.QTc1_median_IQR_Method_1);
x = 1:Gaussian_size;
y = chi2pdf(x,1);

% Normal distribution check: 
[h, p] = chi2gof(QT.QTc1_median_IQR_Method_1,'Alpha',0.05);
[h2, p2] = chi2gof(QT.QTc1_median_IQR_Method_2, 'Alpha',0.05);
[h3, p3] = chi2gof(QT.QTc1_median_human, 'Alpha',0.05);

if h==1 && h2==1 && h3==1
    disp('Chi-squared test: reject null hypothesis, not a normal distribution');
    fprintf('The p-values for the human is %1.3g , Gaussian is %1.3g and the Wavelet is %1.3g .\n',...
        p3, p, p2);
else
    disp('Chi-squared test: null hypotheiss not reject, a normal distribution.');
end

% chi-squared test:
% chi_sqrd = sum((QT.QTc1_median_IQR_Method_1-QT.QTc1_median_IQR_Method_2).^2/...
%     QT.QTc1_median_IQR_Method_1,[1, Gaussian_size]);
% if chi_sqrd == 0
%     disp('No difference between Gaussian and Wavelet')
% else
%     disp('There is a difference by chi-squared method')
% end


%% plot

figure(1)
scatter(x, QT.QTc1_median_IQR_Method_1, 'r', 'filled');
ylabel('Second')
xlabel('Record #')
title('The Gaussian Model distribution for Median QTlc IQR')

figure(2)
x2 = size(QT.QTc1_median_IQR_Method_2);
scatter(1:x2, QT.QTc1_median_IQR_Method_2, 'b', 'filled');
ylabel('Second')
xlabel('Record #')
title('The Wavelet method distribution for Median QTlc IQR')

figure(3)
x3 = size(QT.QTc1_median_human);
scatter(1:x3, QT.QTc1_median_human, 'g', 'filled');
ylabel('Second')
xlabel('Record #')
title('The manual humman annotations for Median QTlc')

% Plot channel 1 and 2 separately

figure(4)
boxplot([QTc1_median_human_ch1 , QTc1_median_method_1_ch1, QTc1_median_method_2_ch1],...
    'Labels',{'Median QTlc Human', 'Median QTlc Gaussian', 'Median QTlc Wavelet'});
ylabel('time (seconds)');
title('QT database Median QTlc for Lead 1')

figure(5)
boxplot([QTc1_median_human_ch2 , QTc1_median_method_1_ch2, QTc1_median_method_2_ch2],...
    'Labels',{'Median QTlc Human', 'Median QTlc Gaussian', 'Median QTlc Wavelet'});
ylabel('time (seconds)');
title('QT database Median QTlc for Lead 2')


% figure(3)
% plot(x,y);
% xlabel('observations')
% ylabel('probability density')
%% Save results for the whole dataset in a CSV


QT_table_colNames = {'subject','record', 'channelNum', 'RR_mean_method_1',...
    'RR_median_method_1', 'QTc1_median_method_1', 'QTc1_mean_method_1',...
    'QTc1_median_IQR_method_1', 'RR_median_method_2', 'RR_mean_method_2','RR_jQRS_method_2',...
    'QTc1_median_method_2', 'QTc1_median_SQI_method_2', 'QTc1_mean_jQRS_method_2',...
    'QTc1_jQRS_method_2_Gauss', 'QTc1_jQRS_method_2_Gauss_SQI', ...
    'QTc1_median_IQR_Method_2', 'RR_median_human', 'RR_mean_human',...
    'QTc1_median_human', 'QTc1_mean_human'};
QT_table = table(QT.subject, QT.record, QT.channelNum, QT.RR_mean_method1,...
    QT.RR_median_method1, QT.QTc1_median_method1, QT.QTc1_mean_method1,...
    QT.QTc1_median_IQR_Method_1, QT.RR_median_method2, QT.RR_mean_method2, QT.RR_jQRS_method2,...
    QT.QTc1_median_method2, QT.QTc1_median_method2_SQI, QT.QTc1_mean_method2_jQRS,...
    QT.QTc1_jQRS_method2_Gauss, QT.QTc1_jQRS_method2_Gauss_SQI,...
    QT.QTc1_median_IQR_Method_2, QT.RR_median_human, QT.RR_mean_human,...
    QT.QTc1_median_human, QT.QTc1_mean_human, 'VariableNames',QT_table_colNames);


csv_fileName = 'QT_wholedataset_results.csv';
fileName = fullfile(results_path, csv_fileName);
writetable(QT_table, fileName);
end