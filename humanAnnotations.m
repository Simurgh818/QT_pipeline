function [humanQT] = humanAnnotations(inPath,annotationFileExtention, fName, figPath, fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Inititalize variables
humanQT.ann={}; humanQT.Tend={}; humanQT.Qstart={}; humanQT.QT={};
humanQT.RR_mean={}; humanQT.RR_median={}; humanQT.QTc1_median={};
humanQT.QTc1_mean={}; humanQT.QTc1_median_IQR={};

% Reading the annotation file
[humanQT.ann{1}, humanQT.ann{2},~,~,humanQT.ann{3}]=rdann(inPath,annotationFileExtention);

% Finding Tend and Qstart
humanQT.Tend = humanQT.ann{1}(humanQT.ann{3}==2);
humanQT.Qstart= humanQT.ann{1}([humanQT.ann{2}=='(' & humanQT.ann{3}==1]);
humanQT.QT = (humanQT.Tend - humanQT.Qstart)/fs;

% ToDo: corrected QT, mean, median of RR and QTc1, and IQR

% making a table for human QT
HumanQT_colNames = {'QT'};

HumanQT_table = table(humanQT.QT, 'VariableNames', HumanQT_colNames);

% save into a csv
csv_fileName = [fName, '_humanQT.csv'];
fileName = fullfile(figPath, csv_fileName);
writetable(HumanQT_table, fileName);

end