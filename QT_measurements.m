function [Lin, Fattahi] = QT_measurements(processedPath, fName, figPath, nChannels, fs)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PQT Interval Estimation:
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

% clear;
close all;

data_base_cor_csv = csvread(processedPath);
fs0 = num2str(fs);
% nChannels_str = num2str(nChannels);
% fs0 = '1000';
% fs = str2double(fs0);
% ToDO: specifiy QT_output.csv location folder
% QT_output = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\QT_pipeline\QT_output.csv';
% [QT_f, RR] = QT_analysis('preprocessed.csv', fs0, nChannels_str, nChannels_str, QT_output, 'w');

% QT = csvread('QT_output.csv', 0,1);
% QT_reshaped = reshape(QT(1:70), [5, 14]);
% md_QT_Qiao = zeros(1,14);

% Initialize variables
QT = zeros(nChannels,5);
RR = zeros(nChannels, 3);
Lin.MeanQT_jQRS=[]; Lin.MedianQT_wavelet=[]; 
Lin.MedianQT_wavelet_SQI=[]; Lin.GaussQT_jQRS=[]; 
Lin.GaussQT_wavelet=[];
Lin.RR_jQRS=[]; Lin.MedianRR_wavelet=[]; Lin.MeanRR_wavelet=[];

Fattahi.MeanQT=[]; Fattahi.MedianQT = []; 
Fattahi.MeanRR=[]; Fattahi.MedianRR = [];
%% 

for ch=1:nChannels
    [QT(ch,:), RR(ch,:)] = QT_analysis_single_lead(data_base_cor_csv(:,ch),fs);

end

%  QT =[meanQT_jQRS, medianQT_wavelet, medianQT_SQI, gaussQT_jQRS, gaussQT_wavelet]

Lin.MeanQT_jQRS = QT(:,1)/fs;
Lin.MedianQT_wavelet = QT(:,2)/fs;
Lin.MedianQT_wavelet_SQI = QT(:,3)/fs;
Lin.GaussQT_jQRS = QT(:,4)/fs;
Lin.GaussQT_wavelet = QT(:,5)/fs;
% QT_cell = {meanQT_jQRS; medianQT_wavelet; medianQT_wavelet_SQI; gaussQT_jQRS...
%     gaussQT_wavelet};

md_QT_Qiao = Lin.MedianQT_wavelet';

% RR = [rr_jQRS, medianRR_wavelet, meanRR_wavelet]
Lin.RR_jQRS = RR(:,1); 
Lin.MedianRR_wavelet =  RR(:,2);
Lin.MeanRR_wavelet =  RR(:,3); 
% RR_cell = {rr_jQRS; medianRR_wavelet, meanRR_wavelet};
%% 2- Model Based method: Mr. Fattahi
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval

% [GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv, fs);
md_QT_Fattahi = zeros(1,nChannels);
rPeaks = [];
qtInt = [];

for ch=1:nChannels
    [GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv(:, ch), fs);
    md_QT_Fattahi(ch) = nanmedian(qtInt(1,:));
    Fattahi.MedianQT(ch,1)= nanmedian(qtInt(1,:));
    Fattahi.MeanQT(ch,1) = nanmean(qtInt(1,:));
    Fattahi.MeanRR(ch,1) = mean(diff(rPeaks)/fs);
    Fattahi.MedianRR(ch,1) = nanmedian(diff(rPeaks))/fs;
end





%% Plotting Dr. Li vs. Mr. Fattahi QT measurments

x = 1:nChannels;
figure(1)
hold on
scatter(x, md_QT_Qiao, 'b', 'filled');
scatter(x, md_QT_Fattahi, 'r', 'filled');
legend('Non-model','Model');
xlabel('channels');
ylabel('time (s)');
hold off
if ~exist(figPath, 'dir')
    mkdir(figPath)
end

figExt = ['.fig';'.eps'; '.png'];
for fx=1:3
    fileName = [fName , figExt(fx,:)];
    figName = fullfile(figPath, fileName);
    saveas(gcf, figName);
end
% close all;

%% save results as a csv

colName = 'channelNum';
Lin_table = struct2table(Lin);
Lin_table.(colName)=  [1:nChannels]';
% Lin_table = addvars(Lin_table, channelNum, 'before', 'MeanQT_jQRS'); 
Fattahi_table = struct2table(Fattahi);
% Fattahi_table = addvars(Fattahi_table, channelNum, 'before', 'MedianQT');
Fattahi_table.(colName) =  [1:nChannels]';
joined_tables = join(Lin_table, Fattahi_table, 'keys', colName);

csv_fileName = [fName, '_QT_results.csv'];
fileName = fullfile(figPath, csv_fileName);
writetable(joined_tables, fileName);
% figure(2)
% x2 = 1:165;
% plot(x2, qtInt(1,1:165,1), x2, qtInt(1,1:165,2))
end
