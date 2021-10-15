% Preprocessing ECG signal
% Sina Dabiri, 2021
clear;
close all;

% Add OSET's Tools folder to the Matlab path
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\OSET\Tools')
% Add Dr. Fattahi's QT folder to Matlab path
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\OSET\UnderDevelopment\QTinterval');
% Add Dr. Qiao's QT estimator folder to Matlab path
%  https://github.com/cliffordlab/QTestimation.git
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\QTestimation\QTestimation\QT_for_Alivecor\')


% Import the PTB dataset using wfdb rdsamp function
[sig, fs, tm] = rdsamp('ptbdb/patient001/s0014lre', 1);
data = sig'; %Transpose the signal to have as row vector.
% plot(tm, sig);
%% Dr. Sameni's OSET baseline removal
% Using the median method to remove baseline:

% figure()
% hold on
% plot(tm, sig, 'b-');

wlen = round(0.6*fs);%set window length as 600 msec
b_md = BaseLine1(data, wlen,'md');

% plot(tm, b_md, 'r-', 'LineWidth',2);
% title("The sliding Median filter");
% hold off

% Using the mean method to remove the baseline:

% figure()
% hold on
% plot(tm, sig, 'b-');
% 'ma' with wlen=0.3 sec
wlen = round(0.8*fs);
b_mn = BaseLine1(data, wlen,'mn');

% plot(tm, b_mn, 'g-', 'LineWidth',2);
% title("The sliding Mean filter");
% hold off

% figure()
data_base_cor = data - b_md;
% plot(tm, data_base_cor, 'b-');
% title("The channel base-corrected");

%% Check to see if Powerline noise is present

% TODO: use periodogram to see if there is a peak at 50 or 60 Hz.

%% Powerline noise cancellation
% notch filter design

fs; % sampling frequency (Hz)
f0 = 60; % powerline frequency (50Hz or 60Hz, depending on the country that
         % the data has been recorded in)
Q = 40; % Q-factor (quality factor) of the notch filter. a parameter that 
        % controls the notch quality and is a compromise between transient 
        % response time and notch quality
Wo = f0/(fs/2); % powerline frequency normalized by the Nyquist frequency 
BW = Wo/Q; % normalized notch filter bandwidth
[b,a] = iirnotch(Wo,BW); % design the filter numerator and denominator
 
% looking at the frequency response
% fvtool(b, a)
x=  data_base_cor;
% filtering a 1-D signal x (apply it on each channel if the signal is 
% multi-dimensional)
y = filter(b, a, x);
% figure()
% plot(tm, y, 'b-')
% title("The channel signal after powerline noise filter");

%% QT Interval Estimation:

% 1- Non-model method: Dr. Qiao wavelet method
%  https://github.com/cliffordlab/QTestimation.git

fs0 = num2str(fs);
y_T = y';
% QT_output = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\Dr. Clifford and Dr. Sameni\QT_pipeline\QT_output.csv'
% QT_analysis(y_T, fs0, '1', '1', 'QT_output.csv', 'a');
y_T = y';
[QT1, RR1] = QT_analysis_single_lead(y_T(:,1),fs) 
md_QT_Qiao = median(QT1)

% 2- Model Based method: Dr. Fattahi
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval

GaussParams=qtParamsGausFit(y_T, fs);
% md_QT_Fattahi = median(GaussParams.q)

% TODO: looks like need to have a multichannel input.
% - what does Fattahi q mean?