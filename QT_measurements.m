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

Troubleshooting = 1; % Troubleshooting flag to plot every step

data_base_cor_csv = csvread(processedPath);%Read preprocessed csv

Method_1.Qon = []; Method_1.Toff=[];
Method_1.QT_mean=[]; Method_1.QT_median = []; 
Method_1.RR_mean=[]; Method_1.RR_median = []; Method_1.QT_median_IQR = [];
Method_1.QTc1_mean=[]; Method_1.QTc1_median=[]; Method_1.QTc1_median_IQR=[];% QTc1 is correction by Sagie's formula
Method_1.QTc2_mean=[]; Method_1.QTc2_median=[]; Method_1.QTc2_median_IQR=[]; % QTc2 is correction by Bazett's formula

QT = zeros(nChannels,5);
RR = zeros(nChannels, 3);
Method_2.Qon=[]; Method_2.Toff=[]; Method_2.Rpeak=[];
Method_2.QR_mean=[]; Method_2.QR_median=[];
Method_2.RT_mean=[]; Method_2.RT_median=[];
Method_2.QT_mean_jQRS=[]; Method_2.QT_median_wavelet=[]; 
Method_2.QT_median_wavelet_SQI=[]; Method_2.GaussQT_jQRS=[]; 
Method_2.GaussQT_SQI=[];
Method_2.RR_jQRS=[]; Method_2.RR_median_wavelet=[]; Method_2.RR_mean_wavelet=[];
Method_2.QT_median_IQR=[]; Method_2.QTc1_mean_jQRS=[]; Method_2.QTc1_median_wavelet=[];
Method_2.QTc1_median_wavelet_SQI=[]; Method_2.GaussQTlc_jQRS=[]; 
Method_2.GaussQTlc_SQI=[]; Method_2.QTc1_median_IQR=[];
%% Dr. Method_2's QT measurment

% Calculating the Q start and T end of the wavelet method
beats = struct();
heasig.nsig=1;

heasig.spf = [1,1];
heasig.spf_ecg = 1;

for ch=1:nChannels
    [QT(ch,:), RR(ch,:), beats] = QT_analysis_single_lead(data_base_cor_csv(:,ch),fs );
%     output R peak from the base level and pass it up every level
%     Calculating Q start and T end
    % resample to 1000Hz
    ecg=resample(data_base_cor_csv(:,ch),1000,fs);
    heasig.freq=1000;
    ecg=lp_filter_1000(ecg);
    heasig.nsamp=length(ecg);
%     beats=wavedet_3D(ecg,[],heasig);
    beats.QRSon = beats.QRSon(~isnan(beats.QRSon));
    beats.Toff = beats.Toff(~isnan(beats.Toff));
    Method_2.Qon(ch,:)= beats.QRSon;
    Method_2.Toff(ch,:)= beats.Toff;
end

% Calculating Q Start and T end
% QR = beats.QRSon;
% RT = beats.Toff;
% Method_2.QR_mean = mean(QR);
% Method_2.RT_mean = mean(RT);


% Since the QT_analysis_single_lead resamples the signal to 1000 Hz, need
% to use this sampling frequency to convert to seconds
fs2 = 1000;
Method_2.QT_mean_jQRS = QT(:,1)/fs2;
Method_2.QT_median_wavelet = QT(:,2)/fs2;
Method_2.QT_median_wavelet_SQI = QT(:,3)/fs2; % SQI stands for signal quality >0.9
% The GaussQT uses a gaussian fit for T wave offset. The jQRS stands for Joachim Behar QRS detector.
% based on Behar Joachim, Jonhson Alistair, Clifford Gari D., Oster Julien A
% Comparison of Single Channel Foetal ECG Extraction Methods. 
% Annals of Biomedical Engineering. 42(6), 1340-53. 2014
Method_2.GaussQT_jQRS = QT(:,4)/fs2; 
Method_2.GaussQT_SQI = QT(:,5)/fs2;



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
    [GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv(:, ch), fs );
%     Quality Control check to exclude QT <0.3 and >0.5 second
    qc_indecies = 0.3<qtInt & qtInt<0.5;
    qtInt = qtInt(qc_indecies);

    Method_1.QT_median(ch,1)= nanmedian(qtInt(1,:));
    Method_1.QT_mean(ch,1) = nanmean(qtInt(1,:));
    Method_1.RR_mean(ch,1) = mean(diff(rPeaks)/fs );
    Method_1.RR_median(ch,1) = nanmedian(diff(rPeaks))/fs ;
end

% Calculate median of interquartile range of 25-75%
Method_1_IQR = iqr(Method_1.QT_median,1);
Method_1_iqr_lower = prctile(Method_1.QT_median, 25)-Method_1_IQR; %calc lower IQR
Method_1_iqr_upper = prctile(Method_1.QT_median, 75)+Method_1_IQR; %calc upper IQR
Method_1_indecies = Method_1_iqr_lower<Method_1.QT_median & Method_1.QT_median<Method_1_iqr_upper;
Method_1.QT_median_IQR(1:nChannels, 1) = median(Method_1.QT_median(Method_1_indecies)); % median QT across the interquartile range

% Correcting QT based on Sagie's Linear regression method: QTlc = QT + 0.154(1-RR) 
% Please note the Wavlet's QRS detector is used for QT correctiohn, the
% Method_2.RR_median_wavelet, instead of the Gaussian's R peak detector,
% which is based on max peak amplitude. 
Method_1.QTc1_mean = Method_1.QT_mean + 0.154*(1-Method_2.RR_mean_wavelet);
Method_1.QTc1_median = Method_1.QT_median + 0.154*(1-Method_2.RR_median_wavelet);
Method_1.QTc1_median_IQR = Method_1.QT_median_IQR + 0.154*(1-Method_2.RR_median_wavelet);


md_QT_Method_1 = Method_1.QTc1_median;

% TODO: Calculate Mr. Method_1 QT_median_SQI
% Dr. Li uses Mbsqi_3 function to select beats, and sets RR <0.3 and >>3 to
% NaN.
%% Quality Control

Method_1.QTc1_median(Method_1.QTc1_median>0.5 | Method_1.QTc1_median<0.3)=NaN;
Method_1.QTc1_mean(Method_1.QTc1_mean>0.5 | Method_1.QTc1_mean<0.3)=NaN;
Method_1.QTc1_median_IQR(Method_1.QTc1_median_IQR>0.5 | Method_1.QTc1_median_IQR<0.3)=NaN;

md_QT_Method_1(md_QT_Method_1>0.5 | md_QT_Method_1<0.3) = NaN;

Method_2.QTc1_median_wavelet(Method_2.QTc1_median_wavelet>0.5 | Method_2.QTc1_median_wavelet<0.3)=NaN;
Method_2.QTc1_median_wavelet_SQI(Method_2.QTc1_median_wavelet_SQI >0.5 | Method_2.QTc1_median_wavelet_SQI<0.3)=NaN;
Method_2.QTc1_median_IQR(Method_2.QTc1_median_IQR >0.5 | Method_2.QTc1_median_IQR<0.3)=NaN;
Method_2.QTc1_mean_jQRS(Method_2.QTc1_mean_jQRS >0.5 | Method_2.QTc1_mean_jQRS<0.3)=NaN;
Method_2.GaussQTlc_jQRS(Method_2.GaussQTlc_jQRS >0.5 | Method_2.GaussQTlc_jQRS<0.3)=NaN;
Method_2.GaussQTlc_SQI(Method_2.GaussQTlc_SQI >0.5 | Method_2.GaussQTlc_SQI<0.3)=NaN;

%% Plotting Dr. Method_2 vs. Mr. Method_1 QT measurments

% x = 1:nChannels;
% figure(1)
% hold on
% scatter(x, md_QT_Qiao, 'b', 'filled');
% scatter(x, md_QT_Method_1, 'r', 'filled');
% legend('wavelet','Gaussian Model');
% title('Median QT');
% xlabel('leads');
% ylabel('time (s)');
% hold off
% if ~exist(figPath, 'dir')
%     mkdir(figPath)
% end
% 
figExt = ['.fig';'.eps'; '.png'];
% for fx=1:3
%     fileName = [fName , figExt(fx,:)];
%     figName = fullfile(figPath, fileName);
%     saveas(gcf, figName);
% end
%% Troubleshooting outliers

% 
[GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv(:, 1), fs);
[L, ~] = size(data_base_cor_csv); %length of the signal

Method_1.Qon = round(rPeaks-abs(waveParams.q(1,:)*fs));% waveParams.q have large negative values fpr recprd se;16272
Method_1.Toff = round(rPeaks+(waveParams.t(1,:)*fs));
Method_1.Qon = Method_1.Qon(~isnan(Method_1.Qon));
Method_1.Qon(Method_1.Qon<0)=0;
Method_1.Qon = Method_1.Qon(Method_1.Qon>0);
Method_1.Toff = Method_1.Toff(~isnan(Method_1.Toff));

Method_2.Rpeak = floor(beats.R(~isnan(beats.R))*fs/fs2); 
Method_2.Qon = round(Method_2.Qon(1,:)*fs/fs2);
Method_2.Toff = round(Method_2.Toff(1,:)*fs/fs2);


% close all;
%% save results as a csv

% colName = 'channelNum';
Method_2_colNames = {'channelNum', 'QTc1_median_wavelet', 'QTc1_median_wavelet_SQI',...
    'GaussQTlc_SQI', 'QTc1_median_IQR', 'QTc1_mean_jQRS','GaussQTlc_jQRS',...
    'RR_jQRS', 'RR_median_wavelet', 'RR_mean_wavelet'};
channelNum = [1:nChannels]';
Method_2_table = table(channelNum, Method_2.QTc1_median_wavelet, Method_2.QTc1_median_wavelet_SQI,...
    Method_2.GaussQTlc_SQI, Method_2.QTc1_median_IQR, Method_2.QTc1_mean_jQRS,... 
    Method_2.GaussQTlc_jQRS, Method_2.RR_jQRS, Method_2.RR_median_wavelet,...
    Method_2.RR_mean_wavelet, 'VariableNames',Method_2_colNames);

Method_1_colNames = {'channelNum', 'QTc1_median', 'QTc1_mean', 'QTc1_median_IQR',...
    'RR_median', 'RR_mean'};

Method_1_table = table(channelNum, Method_1.QTc1_median, Method_1.QTc1_mean,...
    Method_1.QTc1_median_IQR, Method_1.RR_median, Method_1.RR_mean,...
    'VariableNames', Method_1_colNames);

joined_tables = join(Method_1_table, Method_2_table, 'keys', 'channelNum');

csv_fileName = [fName, '_QT_results.csv'];
fileName = fullfile(figPath, csv_fileName);
    writetable(joined_tables, fileName);



% Need to multiply by fs2 
% plot(Method_2.Qon, data_base_cor_csv(Method_2.Qon,1), 'b*', 'LineWidth', 3);
% plot(Method_2.Toff, data_base_cor_csv(Method_2.Toff,1), 'r*', 'LineWidth', 3);
% legend('Lead1','Method 1 Q start', 'Method 1 T end', ...
%     'Method 2 Q start', 'Method 2 T end');
% xlabel('samples');
% ylabel('mV');
% hold off

if Troubleshooting
    figure(3)
    hold on
    plot(1:L, data_base_cor_csv(:,1), 'k-')
    plot(Method_2.Rpeak, data_base_cor_csv(Method_2.Rpeak,1), 'g+')
    xlabel('samples')
    ylabel('mV')
    title("R peak detection step in lead 1")
    legend('preprocessed ', 'Method 2 R peaks')
    hold off

%     figure(4)
%     plot(Method_1.Qon, data_base_cor_csv(Method_1.Qon,1), 'bo', 'LineWidth', 3);
%     plot(Method_1.Toff, data_base_cor_csv(Method_1.Toff,1), 'ro', 'LineWidth', 3);

end
for fx=1:3
    fileName = [fName , '_lead1', figExt(fx,:)];
    figName = fullfile(figPath, fileName);
    saveas(gcf, figName);
end

end
