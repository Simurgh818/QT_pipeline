function QT_measurements(processedPath, fName, figPath, nChannels, fs)
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
% 1- Non-model method: Dr. Qiao wavelet method
%  https://github.com/cliffordlab/QTestimation.git
% Add Dr. Qiao's QT estimator folder to Matlab path
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\QTestimation\QTestimation\QT_for_Alivecor\')
% Add Dr. Fattahi's QT folder to Matlab path
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval
addpath('C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\OSET\UnderDevelopment\QTinterval');

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
nChannels_str = num2str(nChannels);
% fs0 = '1000';
% fs = str2double(fs0);
% ToDO: specifiy QT_output.csv location folder
QT_output = 'C:\Users\sinad\OneDrive - Georgia Institute of Technology\CliffordandSameni\GitRepos\QT_pipeline\QT_output.csv';
% [QT_f, RR] = QT_analysis('preprocessed.csv', fs0, nChannels_str, nChannels_str, QT_output, 'w');

% QT = csvread('QT_output.csv', 0,1);
% QT_reshaped = reshape(QT(1:70), [5, 14]);
% md_QT_Qiao = zeros(1,14);
QT = zeros(nChannels,5);
RR = zeros(nChannels, 3);

for ch=1:nChannels
    [QT(ch,:), RR(ch,:)] = QT_analysis_single_lead(data_base_cor_csv(:,ch),fs);
    
%     md_QT_Qiao(ch) = QT_reshaped(2,ch)/fs; % the second row is 
    % the median of QT values and normalizing by sampling frequency
end
% ToDO: need to normalize by fs 
% restore values in more transparent variables: QT =[meanQT_jQRS,
% medianQT_wavelet, medianQT_SQI, gaussQT_jQRS, gaussQT_wavelet]
% RR = [rr_jQRS, medianRR_wavelet, meanRR_wavelet]

% [QT1, RR1] = QT_analysis_single_lead(data_base_cor_csv(:,1),fs) 
md_QT_Qiao = QT(:,2)'/fs;
%% 2- Model Based method: Mr. Fattahi
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval

[GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(data_base_cor_csv, fs);
md_QT_Fattahi = zeros(1,14);
for ch=1:14
    md_QT_Fattahi(ch) = nanmedian(qtInt(1,:,ch));
end

%% Plotting Dr. Li vs. Mr. Fattahi QT measurments

x = 1:14;
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
close all;

%% save results as a csv

% figure(2)
% x2 = 1:165;
% plot(x2, qtInt(1,1:165,1), x2, qtInt(1,1:165,2))
end
