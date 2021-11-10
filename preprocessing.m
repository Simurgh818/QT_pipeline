function [fs] = preprocessing(inPath, outPath, nChannels)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing ECG signal
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
%% Reading the ECG record using WFDB rdsamp function.

close all;


% Import the PTB dataset using wfdb rdsamp function
[sig, fs, tm] = rdsamp(inPath, 1:nChannels);% for specific channel 

data = sig'; %Transpose the signal to have as row vector.

%% Dr. Sameni's OSET baseline removal
% Using the median method to remove baseline:

wlen = round(0.6*fs);%set window length as 600 msec
b_md = BaseLine1(data, wlen,'md');

% Using the mean method to remove the baseline:

wlen = round(0.8*fs);
b_mn = BaseLine1(data, wlen,'mn');


data_base_cor = data - b_md;


%% Check to see if Powerline noise is present

% TODO: use periodogram to see if there is a peak at 50 or 60 Hz.
% add an if statement to remove power noise if there is a peak
% [L, ~] = size(data_base_cor); %length of the signal
% 
% data_base_cor_f_two_sided = fft(data_base_cor);
% dbc_two_sided = abs(data_base_cor_f_two_sided/L);
% dbc_one_side = dbc_two_sided(1:L/2+1);
% dbc_one_side(2:end-1) = 2*dbc_one_side(2:end-1);
% f = fs*(0:(L/2)) /L;

% figure()
% % periodogram(data_base_cor)]
% plot(f, dbc_one_side)
% title('Single sided amplitute spectrum of base corrected signal')
% xlabel('f(Hz)')
% ylabel('|X|');
% Don't see a peak at 50 or 60 Hz, don't need to do power line noise
% cancellation

%% Powerline noise cancellation
% notch filter design

% fs; % sampling frequency (Hz)
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
% % figure()
% % plot(tm, y, 'b-')
% % title("The channel signal after powerline noise filter");
% figure()
% periodogram(y)
%% CSV write

csvwrite(outPath, y');
end 
