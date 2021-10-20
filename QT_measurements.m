%% QT Interval Estimation:
% Sina Dabiri, 2021

% clear;
close all;
% 1- Non-model method: Dr. Qiao wavelet method
%  https://github.com/cliffordlab/QTestimation.git
% Add Dr. Qiao's QT estimator folder to Matlab path
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\QTestimation\QTestimation\QT_for_Alivecor\')
% Add Dr. Fattahi's QT folder to Matlab path
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\OSET\UnderDevelopment\QTinterval');

data_base_cor_csv = csvread('preprocessed.csv');
fs0 = '1000';
% ToDO: specifiy QT_output.csv location folder
% QT_output = '\QT_pipeline\QT_output.csv'
QT_analysis('preprocessed.csv', fs0, '14', '14', 'QT_output.csv', 'w');

QT = csvread('QT_output.csv', 0,1);
QT_reshaped = reshape(QT(1:70), [5, 14]);
md_QT_Qiao = zeros(1,14);
channel = 1;
for ch=1:14
    md_QT_Qiao(ch) = median(QT_reshaped(:,ch))/fs; % finding the mean of QT values and 
                                         % normalizing by sampling frequency
    
end
% 
% [QT1, RR1] = QT_analysis_single_lead(data_base_cor_csv(:,1),fs) 
% md_QT_Qiao = median(QT1)
%% 2- Model Based method: Dr. Fattahi
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval
fs = str2double(fs0);
[GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv, fs);
md_QT_Fattahi = zeros(1,14);
for ch=1:14
    md_QT_Fattahi(ch) = median(qtInt(1,1:165,ch));
end
md_QT_Fattahi

%% Plotting Dr. Li vs. Mr. Fattahi QT measurments

x = 1:14;
figure(1)
plot(x, md_QT_Qiao, x, md_QT_Fattahi);
legend('Non-model','Model');
xlabel('channels');
ylabel('time (s)');
