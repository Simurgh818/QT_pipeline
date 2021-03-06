function [Method_2, Method_1] = QT_measurements(processedPath, fName, figPath, nChannels, fs, humanQT)
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

Method_1.Qon = {}; Method_1.Toff={}; Method_1.Rpeak={};Method_1.RR={};
Method_1.QT={}; Method_1.QT_mean=[]; Method_1.QT_median = []; 
Method_1.RR_mean=[]; Method_1.RR_median = []; Method_1.QT_median_IQR = [];
Method_1.QTc1_mean=[]; Method_1.QTc1_median=[]; Method_1.QTc1_median_IQR=[];% QTc1 is correction by Sagie's formula
Method_1.QTc2_mean=[]; Method_1.QTc2_median=[]; Method_1.QTc2_median_IQR=[]; % QTc2 is correction by Bazett's formula

QT = zeros(nChannels,5);
RR = zeros(nChannels, 3);
Method_2.Qon={}; Method_2.Toff={}; Method_2.Rpeak={}; Method_2.RR={};
Method_2.RR_mean=[]; Method_2.RR_median=[]; 
Method_2.QT_median=[]; Method_2.QT_mean=[];Method_2.QT={};
Method_2.QT_mean_jQRS=[]; Method_2.QT_median_wavelet=[]; 
Method_2.QT_median_wavelet_SQI=[]; Method_2.GaussQT_jQRS=[]; 
Method_2.GaussQT_SQI=[];
Method_2.RR_jQRS=[]; Method_2.RR_median_wavelet=[]; Method_2.RR_mean_wavelet=[];
Method_2.QTc1_median=[];Method_2.QTc1_mean =[];
Method_2.QT_median_IQR=[]; Method_2.QTc1_mean_jQRS=[]; Method_2.QTc1_median_wavelet=[];
Method_2.QTc1_median_wavelet_SQI=[]; Method_2.GaussQTlc_jQRS=[]; 
Method_2.GaussQTlc_SQI=[]; Method_2.QTc1_median_IQR=[];
%% Dr. Method_2's QT measurment

fs2 = 1000;
[L, ~] = size(data_base_cor_csv); %length of the signal
ecg = [];
ecg_filtered = [];
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
    ecg(:,ch)=resample(data_base_cor_csv(:,ch),fs2,fs);
    resample_back = resample(ecg,fs, fs2);
    ecg_filtered(:,ch)=lp_filter_1000(ecg(:,ch)); % They are using a FIR 
    % linear filter, which has a constat group delay: delay=
    % (length(b)-1)/2, length(b)=233
    D = (233-1)/2;
%     [Xa, Ya, D]=alignsignals(ecg_filtered(:,1), ecg(:,1));
    calc_delay = abs(D*fs/fs2)-1; % manual observation showed a sample difference of 28 instead of 29.
%     Todo: similar to method 1, for each beat that is missing Qon or Toff,
%     sikip the beat
    Method_2.Rpeak(ch,1) = {beats.R};
    Method_2.Qon(ch,1)= {beats.QRSon};
    Method_2.Toff(ch,1)= {beats.Toff};
%     Method_2.Rpeak(ch,:) = {beats.R(~isnan(beats.R))};
    [~, Len_Rpeak] = size(Method_2.Rpeak{ch,1});
    for idx=1 : Len_Rpeak
        if isnan(Method_2.Qon{ch, 1}(idx)) || isnan(Method_2.Toff{ch, 1}(idx)) ||...
                Method_2.Qon{ch, 1}(idx)<0 || Method_2.Toff{ch, 1}(idx)<0
            Method_2.Rpeak{ch, 1}(idx)=NaN ;
            Method_2.Qon{ch, 1}(idx)=NaN;
            Method_2.Toff{ch, 1}(idx)=NaN;
        end
        if  idx >2 && (idx+1)<=Len_Rpeak &&...
                (Method_2.Qon{ch, 1}(idx) >  Method_2.Qon{ch, 1}(idx+1) ||...
                Method_2.Qon{ch, 1}(idx) <  Method_2.Qon{ch, 1}(idx-1))
            Method_2.Qon{ch, 1}(idx)=NaN;
            Method_2.Rpeak{ch, 1}(idx)=NaN;
            Method_2.Toff{ch, 1}(idx)=NaN;
        end

    end
    Method_2.Rpeak{ch,1} = Method_2.Rpeak{ch,1}(~isnan(cell2mat(Method_2.Rpeak(ch,1))));
    Method_2.Qon{ch,1} = Method_2.Qon{ch,1}(~isnan(cell2mat(Method_2.Qon(ch,1))));
    Method_2.Toff{ch,1} = Method_2.Toff{ch,1}(~isnan(cell2mat(Method_2.Toff(ch,1))));

    Method_2.Rpeak(ch,:) = {round(Method_2.Rpeak{ch,:}*fs/fs2-calc_delay)};
    Method_2.RR(ch,:) = {diff(cell2mat(Method_2.Rpeak(ch,:)))/fs};
    Method_2.RR_median(ch,:) = median(cell2mat(Method_2.RR(ch, 1)), 2);
    Method_2.RR_mean(ch,:) = mean(cell2mat(Method_2.RR(ch, 1)), 2);
    
    Method_2.Qon(ch,:) = {round(Method_2.Qon{ch,:}*fs/fs2-calc_delay)};
    Method_2.Toff(ch,:) = {round(Method_2.Toff{ch,:}*fs/fs2-calc_delay)};
    
    Method_2.QT(ch,:) = {(cell2mat(Method_2.Toff(ch,:))-cell2mat(Method_2.Qon(ch,:)))/fs};
    Method_2.QT_median(ch,:) = median(cell2mat(Method_2.QT(ch, :)), 2);
    Method_2.QT_mean(ch,:) = mean(cell2mat(Method_2.QT(ch,:)), 2);
end

if Troubleshooting
%     [L_ecg,~]= size(ecg);
%     figure(1)
%     plot(1:L_ecg,ecg(:,1), 1:L_ecg,ecg_filtered(:,1));
%     legend('resampled ecg', 'filtered ecg');
%     
%     figure(2)
%     plot(1:L_ecg,Xa, 1:L_ecg,Ya(1:L_ecg));
%     legend('resampled ecg', 'aligned filtered ecg');
    
end


% Since the QT_analysis_single_lead resamples the signal to 1000 Hz, need
% to use this sampling frequency to convert to seconds

% ToDo: comment out jQRS and wavelet ones out.
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
Method_2_IQR = iqr(Method_2.QT_median,1);
Method_2_iqr_lower = prctile(Method_2.QT_median, 25)-Method_2_IQR; %calc lower IQR
Method_2_iqr_upper = prctile(Method_2.QT_median, 75)+Method_2_IQR; %calc upper IQR
Method_2_indecies = Method_2_iqr_lower<Method_2.QT_median & Method_2.QT_median <Method_2_iqr_upper;
Method_2.QT_median_IQR(1:nChannels, 1) = median(Method_2.QT_median(Method_2_indecies)); % median QT across the interquartile range

% ToDo: comment out jQRS and wavelet ones out.
% RR = [rr_jQRS, medianRR_wavelet, meanRR_wavelet]
Method_2.RR_jQRS = RR(:,1); 
Method_2.RR_median_wavelet =  RR(:,2);
Method_2.RR_mean_wavelet =  RR(:,3); 


% Correcting QT based on Sagie's Liear regression method: QTlc = QT + 0.154(1-RR)
Method_2.QTc1_median= Method_2.QT_median + 0.154*(1-Method_2.RR_median);
Method_2.QTc1_mean = Method_2.QT_mean + 0.154*(1-Method_2.RR_mean);
Method_2.QTc1_median_IQR = Method_2.QT_median_IQR + 0.154*(1-Method_2.RR_median);

Method_2.QTc1_mean_jQRS = Method_2.QT_mean_jQRS + 0.154*(1-Method_2.RR_jQRS);
Method_2.GaussQTlc_jQRS = Method_2.GaussQT_jQRS + 0.154*(1-Method_2.RR_jQRS);
Method_2.QTc1_median_wavelet = Method_2.QT_median_wavelet + 0.154*(1-Method_2.RR_median_wavelet);
Method_2.QTc1_median_wavelet_SQI = Method_2.QT_median_wavelet_SQI + 0.154*(1-Method_2.RR_median_wavelet);
Method_2.GaussQTlc_SQI = Method_2.GaussQT_SQI + + 0.154*(1-Method_2.RR_median_wavelet);

md_QT_Qiao = Method_2.QTc1_median_wavelet'; %for the plot
%% 2- Model Based method: Mr. Method_1
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval

md_QT_Method_1 = zeros(1,nChannels);
rPeaks = [];
qtInt = [];
GaussParams =[];
soi = [];
waveParams = [];
qc_QT_indices={};

for ch=1:nChannels
%     ToDo: resample to fs of 1000 to have it similar to wavelet method,
%     and it might give more acurate estimates.
    [GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv(:, ch), fs );

    Method_1.Rpeak(ch,:) = {rPeaks};
    Method_1.Qon(ch,:) = {round(rPeaks-abs(waveParams.q(1,:))*fs)};% waveParams.q have large negative values for record se;16272
    Method_1.Toff(ch,:) = {round(rPeaks+(waveParams.t(2,:)*fs))};
    [~, Len_Rpeak] = size(Method_1.Rpeak{ch});
    for idx=1 : Len_Rpeak
        if isnan(Method_1.Qon{ch, 1}(idx)) || isnan(Method_1.Toff{ch, 1}(idx)) ||...
                Method_1.Qon{ch, 1}(idx)<0 || Method_1.Toff{ch, 1}(idx)<0
            Method_1.Rpeak{ch, 1}(idx)=NaN ;
            Method_1.Qon{ch, 1}(idx)=NaN;
            Method_1.Toff{ch, 1}(idx)=NaN;
        end
        if  idx >2 && (idx+1)<=Len_Rpeak &&...
                (Method_1.Qon{ch, 1}(idx) >  Method_1.Qon{ch, 1}(idx+1) ||...
                Method_1.Qon{ch, 1}(idx) <  Method_1.Qon{ch, 1}(idx-1))
            Method_1.Qon{ch, 1}(idx)=NaN;
            Method_1.Rpeak{ch, 1}(idx)=NaN;
            Method_1.Toff{ch, 1}(idx)=NaN;
        end

    end
    Method_1.Rpeak{ch,1} = Method_1.Rpeak{ch,1}(~isnan(cell2mat(Method_1.Rpeak(ch,1))));
    Method_1.Qon{ch,1} = Method_1.Qon{ch,1}(~isnan(cell2mat(Method_1.Qon(ch,1))));
    Method_1.Toff{ch,1} = Method_1.Toff{ch,1}(~isnan(cell2mat(Method_1.Toff(ch,1))));
%     Excluding Qon that don't belong to a beat
    [~, Len_Rpeak] = size(Method_1.Rpeak{ch});
    for idx=1 : Len_Rpeak
        if  idx >2 &&...
                (Method_1.Qon{ch, 1}(idx) >  Method_1.Rpeak{ch, 1}(idx) ||...
                Method_1.Qon{ch, 1}(idx) <  Method_1.Rpeak{ch, 1}(idx-1))
            Method_1.Qon{ch, 1}(idx)=NaN;
            Method_1.Rpeak{ch, 1}(idx)=NaN;
            Method_1.Toff{ch, 1}(idx)=NaN;
        end
        
    end
    Method_1.Rpeak{ch,1} = Method_1.Rpeak{ch,1}(~isnan(cell2mat(Method_1.Rpeak(ch,1))));
    Method_1.Qon{ch,1} = Method_1.Qon{ch,1}(~isnan(cell2mat(Method_1.Qon(ch,1))));
    Method_1.Toff{ch,1} = Method_1.Toff{ch,1}(~isnan(cell2mat(Method_1.Toff(ch,1))));

   %     Quality Control check to exclude QT <0.3 and >0.5 second
%     qc_indecies = 0.3<qtInt & qtInt<0.5;
%     qtInt = qtInt(qc_indecies);
% 
%     Method_1.QT_median(ch,1)= nanmedian(qtInt);
%     Method_1.QT_mean(ch,1) = nanmean(qtInt);
    Method_1.QT(ch,:) = {(cell2mat(Method_1.Toff(ch,:))-cell2mat(Method_1.Qon(ch,:)))/fs};
    qc_QT_indices(ch,:) = {0.3<Method_1.QT{ch,:} & Method_1.QT{ch,:}<0.5};
    Method_1.QT(ch,:) = {Method_1.QT{ch,:}(qc_QT_indices{ch,1})};

    Method_1.QT_median(ch,:) = median(cell2mat(Method_1.QT(ch,:)), 2);
    Method_1.QT_mean(ch,:) = mean(cell2mat(Method_1.QT(ch,:)), 2);

    Method_1.RR(ch,:) = {diff(rPeaks)/fs};
    Method_1.RR_mean(ch,1) = mean(diff(rPeaks)/fs);
    Method_1.RR_median(ch,1) = nanmedian(diff(rPeaks)/fs) ;
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

Method_2.QTc1_median(Method_2.QTc1_median>0.5 | Method_2.QTc1_median<0.3)=NaN;
Method_2.QTc1_mean(Method_2.QTc1_mean>0.5 | Method_2.QTc1_mean<0.3)=NaN;
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

% Method_2.Rpeak = Method_2.Rpeak(1,~isnan(Method_2.Rpeak(1,:)));
% Method_2.Rpeak = resample(Method_2.Rpeak, 1,4);
% Method_2.Rpeak = round(Method_2.Rpeak);
% round(Method_2.Rpeak(1,~isnan())*fs/fs2); 



%% save results as a csv

% colName = 'channelNum';
Method_2_colNames = {'channelNum', 'Method2_QTc1_median', 'Method2_QTc1_mean',...
    'Method2_QTc1_median_IQR', 'Method2_RR_median', 'Method2_RR_mean'};
channelNum = [1:nChannels]';
Method_2_table = table(channelNum, Method_2.QT_median ,Method_2.QTc1_mean,...
    Method_2.QTc1_median_IQR, Method_2.RR_median, Method_2.RR_mean...
    , 'VariableNames',Method_2_colNames);

Method_1_colNames = {'channelNum', 'Method1_QTc1_median', ...
    'Method1_QTc1_mean', 'Method1_QTc1_median_IQR', 'Method1_RR_median',...
    'Method1_RR_mean'};

Method_1_table = table(channelNum, Method_1.QTc1_median, Method_1.QTc1_mean,...
    Method_1.QTc1_median_IQR, Method_1.RR_median, Method_1.RR_mean,...
    'VariableNames', Method_1_colNames);

joined_tables = join(Method_1_table, Method_2_table, 'keys', 'channelNum');

csv_fileName = [fName, '_QT_results.csv'];
fileName = fullfile(figPath, csv_fileName);
    writetable(joined_tables, fileName);

if Troubleshooting
%     figure(3)
%     hold on
%     plot(1:L, data_base_cor_csv(:,1), 'k-');
%     plot(cell2mat(Method_1.Rpeak(1,1)), data_base_cor_csv(cell2mat(Method_1.Rpeak(1,1))), 'm+', 'LineWidth', 3);
%     plot(Method_2.Rpeak(1,:), data_base_cor_csv(Method_2.Rpeak(1,:),1), 'g+','LineWidth', 3);
%     plot(humanQT.R, data_base_cor_csv(humanQT.R,1), 'r+','LineWidth', 3);
%     xlabel('samples')
%     ylabel('mV')
%     title("R peak detection step in lead 1")
%     legend({'preprocessed', 'Method 1 R peaks', 'Method 2 R peaks', 'manual R peaks'});
%     hold off
    qc_QT_indices_ch1 = cell2mat(qc_QT_indices(1,1));
    Method_1.Qon_ch1 = Method_1.Qon{1,:}(qc_QT_indices_ch1);
    figure(4)
    hold on
    plot(1:L, data_base_cor_csv(:,1), 'k-');
    plot(Method_1.Qon_ch1 , data_base_cor_csv(Method_1.Qon_ch1 ,1), 'mo', 'LineWidth', 3);
    plot(cell2mat(Method_1.Rpeak(1,1)), data_base_cor_csv(cell2mat(Method_1.Rpeak(1,1))), 'b+', 'LineWidth', 3);
    plot(cell2mat(Method_2.Qon(1,1)), data_base_cor_csv(cell2mat(Method_2.Qon(1,1)),1), 'go', 'LineWidth', 3);
    plot(humanQT.Qstart, data_base_cor_csv(humanQT.Qstart,1), 'ro','LineWidth', 3);
    legend('Preprocessed ECG','Gaussian Qon', 'Gaussian R peaks', 'Wavelet Qon', 'manual Qon')
    title("Qon detection step in lead 1")
    xlabel('time');
    ylabel('[mV]');
    hold off

    figure(5)
    hold on
    plot(1:L, data_base_cor_csv(:,1), 'k-');
    plot(cell2mat(Method_1.Toff(1,1)), data_base_cor_csv(cell2mat(Method_1.Toff(1,1)),1), 'm*', 'LineWidth', 3);
    plot(cell2mat(Method_2.Toff(1,1)), data_base_cor_csv(cell2mat(Method_2.Toff(1,1)),1), 'g*', 'LineWidth', 3);
    plot(humanQT.Tend, data_base_cor_csv(humanQT.Tend,1), 'r*','LineWidth', 3);
    legend('Preprocessed ECG','Gaussian Toff', 'Wavelet Toff', 'manual Toff')
    title("Toff detection step in lead 1")
    xlabel('time');
    ylabel('[mV]');
    hold off

end
% close all;
% for fx=1:3
%     fileName = [fName , '_lead1', figExt(fx,:)];
%     figName = fullfile(figPath, fileName);
%     saveas(gcf, figName);
% end

end
