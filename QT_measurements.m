function [Li, Fattahi] = QT_measurements(processedPath, fName, figPath, nChannels, fs)
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
Li.MeanQT_jQRS=[]; Li.MedianQT_wavelet=[]; 
Li.MedianQT_wavelet_SQI=[]; Li.GaussQT_jQRS=[]; 
Li.GaussQT_wavelet=[];
Li.RR_jQRS=[]; Li.MedianRR_wavelet=[]; Li.MeanRR_wavelet=[];
Li.MedianQT_IQR=[]; Li.MeanQTlc_jQRS=[]; Li.MedianQTlc_wavelet=[];
Li.MedianQTlc_wavelet_SQI=[]; Li.GaussQTlc_jQRS=[]; 
Li.GaussQTlc_wavelet=[]; Li.MedianQTlc_IQR=[];

Fattahi.MeanQT=[]; Fattahi.MedianQT = []; 
Fattahi.MeanRR=[]; Fattahi.MedianRR = []; Fattahi.MedianQT_IQR = [];
Fattahi.MeanQTlc=[]; Fattahi.MedianQTlc=[]; Fattahi.MedianQTlc_IQR=[];
%% Dr. Li's QT measurment

for ch=1:nChannels
    [QT(ch,:), RR(ch,:)] = QT_analysis_single_lead(data_base_cor_csv(:,ch),fs);

end

Li.MeanQT_jQRS = QT(:,1)/fs;
Li.MedianQT_wavelet = QT(:,2)/fs;
Li.MedianQT_wavelet_SQI = QT(:,3)/fs;
Li.GaussQT_jQRS = QT(:,4)/fs;
Li.GaussQT_wavelet = QT(:,5)/fs;

md_QT_Qiao = Li.MedianQT_wavelet'; %for the plot

% Calculate median of interquartile range of 25-75%
Li_IQR = iqr(Li.MedianQT_wavelet,1);
Li_iqr_lower = prctile(Li.MedianQT_wavelet, 25)-Li_IQR; %calc lower IQR
Li_iqr_upper = prctile(Li.MedianQT_wavelet, 75)+Li_IQR; %calc upper IQR
Li_indecies = Li_iqr_lower<Li.MedianQT_wavelet & Li.MedianQT_wavelet <Li_iqr_upper;
Li.MedianQT_IQR(1:nChannels, 1) = median(Li.MedianQT_wavelet(Li_indecies)); % median QT across the interquartile range

% RR = [rr_jQRS, medianRR_wavelet, meanRR_wavelet]
Li.RR_jQRS = RR(:,1); 
Li.MedianRR_wavelet =  RR(:,2);
Li.MeanRR_wavelet =  RR(:,3); 

% Correcting QT based on Sagie's Liear regression method: QTlc = QT + 0.154(1-RR) 
Li.MeanQTlc_jQRS = Li.MeanQT_jQRS + 0.154*(1-Li.RR_jQRS);
Li.GaussQTlc_jQRS = Li.GaussQT_jQRS + 0.154*(1-Li.RR_jQRS);

Li.MedianQTlc_wavelet = Li.MedianQT_wavelet + 0.154*(1-Li.MedianRR_wavelet);
Li.MedianQTlc_wavelet_SQI = Li.MedianQT_wavelet_SQI + 0.154*(1-Li.MedianRR_wavelet);
Li.GaussQTlc_wavelet = Li.GaussQT_wavelet + + 0.154*(1-Li.MedianRR_wavelet);
Li.MedianQTlc_IQR = Li.MedianQT_IQR + 0.154*(1-Li.MedianRR_wavelet);
%% 2- Model Based method: Mr. Fattahi
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval

md_QT_Fattahi = zeros(1,nChannels);
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

% Correcting QT based on Sagie's Linear regression method: QTlc = QT + 0.154(1-RR) 
Fattahi.MeanQTlc = Fattahi.MeanQT + 0.154*(1-Fattahi.MeanRR);
Fattahi.MedianQTlc = Fattahi.MedianQT + 0.154*(1-Fattahi.MedianRR);
Fattahi.MedianQTlc_IQR = Fattahi.MedianQT_IQR + 0.154*(1-Fattahi.MedianRR);

% TODO: Calculate Mr. Fattahi MedianQT_SQI
% Dr. Li uses Mbsqi_3 function to select beats, and sets RR <0.3 and >>3 to
% NaN.
%% Quality Control

Fattahi.MedianQTlc(Fattahi.MedianQTlc>1 | Fattahi.MedianQTlc<0.1)=NaN;
Fattahi.MeanQTlc(Fattahi.MeanQTlc>1 | Fattahi.MeanQTlc<0.1)=NaN;

md_QT_Fattahi(md_QT_Fattahi>1 | md_QT_Fattahi<0.1) = NaN;
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
Li_colNames = {'channelNum', 'MedianQTlc_wavelet', 'MedianQTlc_wavelet_SQI',...
    'GaussQTlc_wavelet', 'MedianQTlc_IQR', 'MeanQTlc_jQRS','GaussQTlc_jQRS',...
    'RR_jQRS', 'MedianRR_wavelet'};
channelNum = [1:nChannels]';
Li_table = table(channelNum, Li.MedianQTlc_wavelet, Li.MedianQTlc_wavelet_SQI,...
    Li.GaussQTlc_wavelet, Li.MedianQTlc_IQR, Li.MeanQTlc_jQRS,... 
    Li.GaussQTlc_jQRS, Li.RR_jQRS, Li.MedianRR_wavelet,'VariableNames',Li_colNames);

Fattahi_colNames = {'channelNum', 'MedianQTlc', 'MeanQTlc', 'MedianQTlc_IQR',...
    'MedianRR', 'MeanRR'};

Fattahi_table = table(channelNum, Fattahi.MedianQTlc, Fattahi.MeanQTlc,...
    Fattahi.MedianQTlc_IQR, Fattahi.MedianRR, Fattahi.MeanRR,...
    'VariableNames', Fattahi_colNames);

joined_tables = join(Fattahi_table, Li_table, 'keys', 'channelNum');

csv_fileName = [fName, '_QT_results.csv'];
fileName = fullfile(figPath, csv_fileName);
writetable(joined_tables, fileName);

end
