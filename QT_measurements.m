function [Method_2, Method_1] = QT_measurements(processedPath, fName, figPath, nChannels, fs)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QT Interval Estimation:
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
%% Initialize variables

data_base_cor_csv = csvread(processedPath);%Read preprocessed csv

QT = zeros(nChannels,5);
RR = zeros(nChannels, 3);
Method_2.QT_mean_jQRS=[]; Method_2.QT_median_wavelet=[]; 
Method_2.QT_median_wavelet_SQI=[]; Method_2.GaussQT_jQRS=[]; 
Method_2.GaussQT_SQI=[];
Method_2.RR_jQRS=[]; Method_2.RR_median_wavelet=[]; Method_2.RR_mean_wavelet=[];
Method_2.QT_median_IQR=[]; Method_2.QTc1_mean_jQRS=[]; Method_2.QTc1_median_wavelet=[];
Method_2.QTc1_median_wavelet_SQI=[]; Method_2.GaussQTlc_jQRS=[]; 
Method_2.GaussQTlc_SQI=[]; Method_2.QTc1_median_IQR=[];

Method_1.QT_mean=[]; Method_1.QT_median = []; 
Method_1.RR_mean=[]; Method_1.RR_median = []; Method_1.QT_median_IQR = [];
Method_1.QTc1_mean=[]; Method_1.QTc1_median=[]; Method_1.QTc1_median_IQR=[];
%% Dr. Method_2's QT measurment

for ch=1:nChannels
    [QT(ch,:), RR(ch,:)] = QT_analysis_single_lead(data_base_cor_csv(:,ch),fs);

end

Method_2.QT_mean_jQRS = QT(:,1)/fs;
Method_2.QT_median_wavelet = QT(:,2)/fs;
Method_2.QT_median_wavelet_SQI = QT(:,3)/fs;
Method_2.GaussQT_jQRS = QT(:,4)/fs;
Method_2.GaussQT_SQI = QT(:,5)/fs;



% Calculate median of interquartile range of 25-75%
Method_2_IQR = iqr(Method_2.QT_median_wavelet,1);
Method_2_iqr_lower = prctile(Method_2.QT_median_wavelet, 25)-Method_2_IQR; %calc lower IQR
Method_2_iqr_upper = prctile(Method_2.QT_median_wavelet, 75)+Method_2_IQR; %calc upper IQR
Method_2_indecies = Method_2_iqr_lower<Method_2.QT_median_wavelet & Method_2.QT_median_wavelet <Method_2_iqr_upper;
Method_2.QT_median_IQR(1:nChannels, 1) = median(Method_2.QT_median_wavelet(Method_2_indecies)); % median QT across the interquartile range

% RR = [rr_jQRS, medianRR_wavelet, meanRR_wavelet]
Method_2.RR_jQRS = RR(:,1); 
Method_2.RR_median_wavelet =  RR(:,2);
Method_2.RR_mean_wavelet =  RR(:,3); 

% Correcting QT based on Sagie's Liear regression method: QTlc = QT + 0.154(1-RR) 
Method_2.QTc1_mean_jQRS = Method_2.QT_mean_jQRS + 0.154*(1-Method_2.RR_jQRS);
Method_2.GaussQTlc_jQRS = Method_2.GaussQT_jQRS + 0.154*(1-Method_2.RR_jQRS);

Method_2.QTc1_median_wavelet = Method_2.QT_median_wavelet + 0.154*(1-Method_2.RR_median_wavelet);
Method_2.QTc1_median_wavelet_SQI = Method_2.QT_median_wavelet_SQI + 0.154*(1-Method_2.RR_median_wavelet);
Method_2.GaussQTlc_SQI = Method_2.GaussQT_SQI + + 0.154*(1-Method_2.RR_median_wavelet);
Method_2.QTc1_median_IQR = Method_2.QT_median_IQR + 0.154*(1-Method_2.RR_median_wavelet);

md_QT_Qiao = Method_2.QTc1_median_wavelet'; %for the plot
%% 2- Model Based method: Mr. Method_1
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval

md_QT_Method_1 = zeros(1,nChannels);
rPeaks = [];
qtInt = [];
GaussParams =[];
soi = [];
waveParams = [];

for ch=1:nChannels
    [GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv(:, ch), fs);
%     Quality Control check to exclude QT <0.1 and >1 second
    qc_indecies = 0.1<qtInt & qtInt<1;
    qtInt = qtInt(qc_indecies);

    Method_1.QT_median(ch,1)= nanmedian(qtInt(1,:));
    Method_1.QT_mean(ch,1) = nanmean(qtInt(1,:));
    Method_1.RR_mean(ch,1) = mean(diff(rPeaks)/fs);
    Method_1.RR_median(ch,1) = nanmedian(diff(rPeaks))/fs;
end

% Calculate median of interquartile range of 25-75%
Method_1_IQR = iqr(Method_1.QT_median,1);
Method_1_iqr_lower = prctile(Method_1.QT_median, 25)-Method_1_IQR; %calc lower IQR
Method_1_iqr_upper = prctile(Method_1.QT_median, 75)+Method_1_IQR; %calc upper IQR
Method_1_indecies = Method_1_iqr_lower<Method_1.QT_median & Method_1.QT_median<Method_1_iqr_upper;
Method_1.QT_median_IQR(1:nChannels, 1) = median(Method_1.QT_median(Method_1_indecies)); % median QT across the interquartile range

% Correcting QT based on Sagie's Linear regression method: QTlc = QT + 0.154(1-RR) 
Method_1.QTc1_mean = Method_1.QT_mean + 0.154*(1-Method_1.RR_mean);
Method_1.QTc1_median = Method_1.QT_median + 0.154*(1-Method_1.RR_median);
Method_1.QTc1_median_IQR = Method_1.QT_median_IQR + 0.154*(1-Method_1.RR_median);


md_QT_Method_1 = Method_1.QTc1_median;

% TODO: Calculate Mr. Method_1 QT_median_SQI
% Dr. Li uses Mbsqi_3 function to select beats, and sets RR <0.3 and >>3 to
% NaN.
%% Quality Control

Method_1.QTc1_median(Method_1.QTc1_median>1 | Method_1.QTc1_median<0.1)=NaN;
Method_1.QTc1_mean(Method_1.QTc1_mean>1 | Method_1.QTc1_mean<0.1)=NaN;

md_QT_Method_1(md_QT_Method_1>1 | md_QT_Method_1<0.1) = NaN;
%% Plotting Dr. Method_2 vs. Mr. Method_1 QT measurments

x = 1:nChannels;
figure(1)
hold on
scatter(x, md_QT_Qiao, 'b', 'filled');
scatter(x, md_QT_Method_1, 'r', 'filled');
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
% close all;
%% save results as a csv

% colName = 'channelNum';
Method_2_colNames = {'channelNum', 'QTc1_median_wavelet', 'QTc1_median_wavelet_SQI',...
    'GaussQTlc_SQI', 'QTc1_median_IQR', 'QTc1_mean_jQRS','GaussQTlc_jQRS',...
    'RR_jQRS', 'RR_median_wavelet'};
channelNum = [1:nChannels]';
Method_2_table = table(channelNum, Method_2.QTc1_median_wavelet, Method_2.QTc1_median_wavelet_SQI,...
    Method_2.GaussQTlc_SQI, Method_2.QTc1_median_IQR, Method_2.QTc1_mean_jQRS,... 
    Method_2.GaussQTlc_jQRS, Method_2.RR_jQRS, Method_2.RR_median_wavelet,'VariableNames',Method_2_colNames);

Method_1_colNames = {'channelNum', 'QTc1_median', 'QTc1_mean', 'QTc1_median_IQR',...
    'RR_median', 'RR_mean'};

Method_1_table = table(channelNum, Method_1.QTc1_median, Method_1.QTc1_mean,...
    Method_1.QTc1_median_IQR, Method_1.RR_median, Method_1.RR_mean,...
    'VariableNames', Method_1_colNames);

joined_tables = join(Method_1_table, Method_2_table, 'keys', 'channelNum');

csv_fileName = [fName, '_QT_results.csv'];
fileName = fullfile(figPath, csv_fileName);
writetable(joined_tables, fileName);

end
