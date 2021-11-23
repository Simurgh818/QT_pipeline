function [humanQT] = humanAnnotations(inPath,annotationFileExtention, fName, figPath, fs, nChannels)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Inititalize variables
humanQT.ann={}; humanQT.Tend={}; humanQT.Qstart={}; humanQT.QT=[];
humanQT.RR=[] ; humanQT.RR_mean=[]; humanQT.RR_median=[]; humanQT.QTc1_median=[];
humanQT.QTc1_mean=[]; 

% Reading the annotation file
[humanQT.ann{1}, humanQT.ann{2},~,~,humanQT.ann{3}]=rdann(inPath,annotationFileExtention);

% Finding Tend and Qstart
humanQT.Tend = humanQT.ann{1}([humanQT.ann{2}==')' & humanQT.ann{3}==2]);
humanQT.Qstart= humanQT.ann{1}([humanQT.ann{2}=='(' & humanQT.ann{3}==1]);
[dimQ,~] = size(humanQT.Qstart);
if isempty(humanQT.Tend) || isempty(humanQT.Qstart)
    humanQT.QT = NaN;
else
    for q=1:dimQ
        humanQT.QT(q,1) = (humanQT.Tend(q) - humanQT.Qstart(q))/fs;
    end
end

% Finding RR interval
humanQT.R=humanQT.ann{1}([humanQT.ann{2}=='N']);
[dimQT, ~] = size(humanQT.R);
for r=1:(dimQT-1)
    humanQT.RR(r,1) = (humanQT.R(r+1)-humanQT.R(r))/fs;
end

% median RR
humanQT.RR_median = median(humanQT.RR);
humanQT.RR_mean = mean(humanQT.RR);

% Correcting QT based on Sagie's Liear regression method: QTlc = QT + 0.154(1-RR) 
humanQT.QTc1  = humanQT.QT + 0.154*(1-humanQT.RR_median);
humanQT.QTc1_median = median(humanQT.QTc1);
humanQT.QTc1_mean = mean(humanQT.QTc1);

% Since for the QT dataset the annotation was done based on both channels
for ch=1:nChannels
    humanQT.RR_median(ch,1) = humanQT.RR_median;
    humanQT.RR_mean(ch,1) = humanQT.RR_mean;
    humanQT.QTc1_median(ch,1) = humanQT.QTc1_median;
    humanQT.QTc1_mean(ch,1) = humanQT.QTc1_mean;
end
%% Quality Control
humanQT.QTc1_median(humanQT.QTc1_median>0.5 | humanQT.QTc1_median<0.3)=NaN;
humanQT.QTc1_mean(humanQT.QTc1_mean>0.5 | humanQT.QTc1_mean<0.3)=NaN;

%% making a table for human QT
channelNum = [1:nChannels]';
HumanQT_colNames = {'channelNum', 'RR_median_human', 'RR_mean_human', 'QTc1_median_human',...
    'QTc1_mean_human'};

HumanQT_table = table(channelNum, humanQT.RR_median, humanQT.RR_mean,...
    humanQT.QTc1_median, humanQT.QTc1_mean , 'VariableNames', HumanQT_colNames);

% save into a csv
csv_fileName = [fName, '_humanQT.csv'];
fileName = fullfile(figPath, csv_fileName);
writetable(HumanQT_table, fileName);

end