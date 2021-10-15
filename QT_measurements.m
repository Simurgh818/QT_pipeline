%% QT Interval Estimation:
% Sina Dabiri, 2021


% 1- Non-model method: Dr. Qiao wavelet method
%  https://github.com/cliffordlab/QTestimation.git

fs0 = num2str(fs);
data_base_cor_T = data_base_cor';
% QT_output = '\QT_pipeline\QT_output.csv'
% QT_analysis(data_base_cor_T, fs0, '1', '1', 'QT_output.csv', 'a');

[QT1, RR1] = QT_analysis_single_lead(data_base_cor_T(:,1),fs) 
md_QT_Qiao = median(QT1)

% 2- Model Based method: Dr. Fattahi
% https://github.com/alphanumericslab/OSET/tree/master/UnderDevelopment/QTinterval

GaussParams=qtParamsGausFit(data_base_cor_T, fs);
% md_QT_Fattahi = median(GaussParams.q)

% TODO: looks like need to have a multichannel input.
% - what does Fattahi q mean?