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
% close all;

data_base_cor_csv = csvread(processedPath);
% fs0 = num2str(fs);
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
Lin.MedianQT_IQR=[]; Lin.MeanQTlc_jQRS=[]; Lin.MedianQTlc_wavelet=[];
Lin.MedianQTlc_wavelet_SQI=[]; Lin.GaussQTlc_jQRS=[]; 
Lin.GaussQTlc_wavelet=[]; Lin.MedianQTlc_IQR=[];

Fattahi.MeanQT=[]; Fattahi.MedianQT = []; 
Fattahi.MeanRR=[]; Fattahi.MedianRR = []; Fattahi.MedianQT_IQR = [];
Fattahi.MeanQTlc=[]; Fattahi.MedianQTlc=[]; Fattahi.MedianQTlc_IQR=[];
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
% Calculate median of interquartile range of 25-75%
Lin_IQR = iqr(Lin.MedianQT_wavelet,1);
Lin_iqr_lower = prctile(Lin.MedianQT_wavelet, 25)-Lin_IQR; %calc lower IQR
Lin_iqr_upper = prctile(Lin.MedianQT_wavelet, 75)+Lin_IQR; %calc upper IQR
Lin_indecies = Lin_iqr_lower<Lin.MedianQT_wavelet & Lin.MedianQT_wavelet <Lin_iqr_upper;
Lin.MedianQT_IQR(1:nChannels, 1) = median(Lin.MedianQT_wavelet(Lin_indecies)); % median QT across the interquartile range

% RR = [rr_jQRS, medianRR_wavelet, meanRR_wavelet]
Lin.RR_jQRS = RR(:,1); 
Lin.MedianRR_wavelet =  RR(:,2);
Lin.MeanRR_wavelet =  RR(:,3); 

% Correcting QT based on Sagie's linear regression method: QTlc = QT + 0.154(1-RR) 
Lin.MeanQTlc_jQRS = Lin.MeanQT_jQRS + 0.154*(1-Lin.RR_jQRS);
Lin.GaussQTlc_jQRS = Lin.GaussQT_jQRS + 0.154*(1-Lin.RR_jQRS);

Lin.MedianQTlc_wavelet = Lin.MedianQT_wavelet + 0.154*(1-Lin.MedianRR_wavelet);
Lin.MedianQTlc_wavelet_SQI = Lin.MedianQT_wavelet_SQI + 0.154*(1-Lin.MedianRR_wavelet);
Lin.GaussQTlc_wavelet = Lin.GaussQT_wavelet + + 0.154*(1-Lin.MedianRR_wavelet);
Lin.MedianQTlc_IQR = Lin.MedianQT_IQR + 0.154*(1-Lin.MedianRR_wavelet);
%% 2- Model Based method: Mr. Fattahi
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval

% [GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv, fs);
md_QT_Fattahi = zeros(1,nChannels);
rPeaks = [];
qtInt = [];
GaussParams =[];
soi = [];
waveParams = [];

for ch=1:nChannels
    [GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv(:, ch), fs);
    md_QT_Fattahi(ch) = nanmedian(qtInt(1,:));
    Fattahi.MedianQT(ch,1)= nanmedian(qtInt(1,:));
    Fattahi.MeanQT(ch,1) = nanmean(qtInt(1,:));
    Fattahi.MeanRR(ch,1) = mean(diff(rPeaks)/fs);
    Fattahi.MedianRR(ch,1) = nanmedian(diff(rPeaks))/fs;
end

% Calculate median of interquartile range of 25-75%
Fattahi_IQR = iqr(Fattahi.MedianQT,1);
Fattahi_iqr_lower = prctile(Fattahi.MedianQT, 25)-Fattahi_IQR; %calc lower IQR
Fattahi_iqr_upper = prctile(Fattahi.MedianQT, 75)+Fattahi_IQR; %calc upper IQR
Fattahi_indecies = Fattahi_iqr_lower<Fattahi.MedianQT & Fattahi.MedianQT<Fattahi_iqr_upper;
Fattahi.MedianQT_IQR(1:nChannels, 1) = median(Fattahi.MedianQT(Fattahi_indecies)); % median QT across the interquartile range

% Correcting QT based on Sagie's linear regression method: QTlc = QT + 0.154(1-RR) 
Fattahi.MeanQTlc = Fattahi.MeanQT + 0.154*(1-Fattahi.MeanRR);
Fattahi.MedianQTlc = Fattahi.MedianQT + 0.154*(1-Fattahi.MedianRR);
Fattahi.MedianQTlc_IQR = Fattahi.MedianQT_IQR + 0.154*(1-Fattahi.MedianRR);

% TODO: Calculate Mr. Fattahi MedianQT_SQI

%% Plotting Dr. Li vs. Mr. Fattahi QT measurments

x = 1:nChannels;
figure(1)
hold on
scatter(x, md_QT_Qiao, 'b', 'filled');
scatter(x, md_QT_Fattahi, 'r', 'filled');
legend('wavelet','Gaussian Model');
title('Median QT');
xlabel('leads');
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
%% Troubleshooting outliers

% 
% [GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv(:, 9), fs);
% [L, ~] = size(data_base_cor_csv); %length of the signal
% 
% Q_start = round(rPeaks+waveParams.q(1,:)*fs);
% T_end = round(rPeaks+waveParams.t(2,:)*fs);
% Q_start = Q_start(Q_start>0);
% T_end = T_end(~isnan(T_end));
% figure(2)
% hold on
% plot(1:L, -data_base_cor_csv(:,9), 'b-');
% title('Lead 9');
% plot(Q_start, -data_base_cor_csv(Q_start,9), 'g*');
% plot(T_end, -data_base_cor_csv(T_end,9), 'r*');
% legend('Lead 9','Q start', 'T end');
% xlabel('samples');
% ylabel('mV');
% hold off
% for fx=1:3
%     fileName = [fName , '_lead9', figExt(fx,:)];
%     figName = fullfile(figPath, fileName);
%     saveas(gcf, figName);
% end
close all;
%% save results as a csv

% colName = 'channelNum';
Lin_colNames = {'channelNum', 'MedianQTlc_wavelet', 'MedianQTlc_wavelet_SQI',...
    'GaussQTlc_wavelet', 'MedianQTlc_IQR', 'MeanQTlc_jQRS','GaussQTlc_jQRS',...
    'RR_jQRS', 'MedianRR_wavelet'};
channelNum = [1:nChannels]';
Lin_table = table(channelNum, Lin.MedianQTlc_wavelet, Lin.MedianQTlc_wavelet_SQI,...
    Lin.GaussQTlc_wavelet, Lin.MedianQTlc_IQR, Lin.MeanQTlc_jQRS,... 
    Lin.GaussQTlc_jQRS, Lin.RR_jQRS, Lin.MedianRR_wavelet,'VariableNames',Lin_colNames);

% Lin_table = struct2table(Lin);
% Lin_table.(colName)=  [1:nChannels]';
% Lin_table = addvars(Lin_table, channelNum, 'before', 'MeanQT_jQRS');
Fattahi_colNames = {'channelNum', 'MedianQTlc', 'MeanQTlc', 'MedianQTlc_IQR',...
    'MedianRR', 'MeanRR'};

Fattahi_table = table(channelNum, Fattahi.MedianQTlc, Fattahi.MeanQTlc,...
    Fattahi.MedianQTlc_IQR, Fattahi.MedianRR, Fattahi.MeanRR,...
    'VariableNames', Fattahi_colNames);
% Fattahi_table = addvars(Fattahi_table, channelNum, 'before', 'MedianQT');
% Fattahi_table.(colName) =  [1:nChannels]';
joined_tables = join(Fattahi_table, Lin_table, 'keys', 'channelNum');

csv_fileName = [fName, '_QT_results.csv'];
fileName = fullfile(figPath, csv_fileName);
writetable(joined_tables, fileName);
% figure(2)
% x2 = 1:165;
% plot(x2, qtInt(1,1:165,1), x2, qtInt(1,1:165,2))
end
