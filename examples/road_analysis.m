%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM INPUT DATA ANALYSIS
% --------------------------
% Descr.: Analysis methods and techniques for 
%         Toyota track experiments with unknown input
% Author: Thomas Beaudiun, Shota Yamada
%         University of Tokyo, HFlab, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc; 
load('data/woctrl_2.mat');
x = aux_rr(1:roundx(length(aux_rr),-2));
load('data/woctrl_4.mat');
x = [x;aux_rr(1:roundx(length(aux_rr),-2))];
fs = floor(1/time2(2));

%% Wavelet
harm.df = 0.1;                  % freq resolution
harm.fl = 1;                    % lowest freq
harm.fh = 60;                   % highest freq
wave.basis = 'morlet';          % wavelet type
wave.figs = 0;                  % output boolean
sigc.level = 0.95;              % significance level
sigc.alpha = 0.1;               % lag coefficient
sigc.range = [5,30];
[power,confd,Fw,time] = pwavelet(x, fs, harm, wave, sigc);
Pwn = (power.f-mean(power.f))./std(power.f);

%% PSD
order = 16;
harm.df = 0.1;                  % freq resolution
wd = fs/harm.df;
[Px,Fx] = pwelch(x,hanning(wd),wd/2,Fw,fs);
Pxn = (dbm(Px)-mean(dbm(Px)))./std(dbm(Px));
[Pl,Fl] = pyulear(x,order,Fw,fs);
Pln = (dbm(Pl)-mean(dbm(Pl)))./std(dbm(Pl));

%% Short-Time Fourier Transform (STFT)
[time,Ft,power,avg] = pstft(time,x,[1,60],0.1,0.2,[5,30]);
Ptn = (avg-mean(avg))./std(avg);

hfig=figure;
semilogx(Fx,Pxn), hold on
semilogx(Fw,Pwn,Fl,Pln,Ft,Ptn)
    xlim([harm.fl,harm.fh]);
    legend('psd','wavelet','welch','stft')
hfig = pubfig(hfig);
    hfig.LineWidth = [2.2,2.2,2.2,1];