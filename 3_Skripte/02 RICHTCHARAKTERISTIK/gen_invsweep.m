function [sinv] = gen_invsweep(T,f1,f2,fs,max_sig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Create the swept sine tone
w1 = 2*pi*f1; %Omega1
w2 = 2*pi*f2; %Omega2
K = T*w1/log(w2/w1); %growingparameter for creating the sweep
L = T/log(w2/w1); %lengthparameter for creating the sweep
t = linspace(0,T-1/fs,fs*T); %time vector
s = sin(K*(exp(t/L) - 1)); %exponential sweep timesignal
% LT=length(s); %length of sweep vector

sinv=s(length(s):-1:1).*exp(-t./L);%generate backwards exponential sine sweep


% normalization occurs in main skript.
%% find max value for filtered signal to normalize invsweep in time dimension.
%% maximum = maxsig;
% sinv=invsweep/max(abs(invsweep)); %normalize timesignal to maximum value of time signal
% sinv=invsweep*maximum; %normalize timesignal to maximum value of time signal



% PLOT SINV
% NFFT=145576
% SINV=fft(sinv)
% f = FS/2*linspace(0,1,NFFT/2); %create plotting vector
% loglog(f,abs(SINV(1:NFFT/2))/max(abs(SINV))) %Plot single-sidedamplitude spectrum of inverse filter
% title('Amplitude Spectrum of the Sine Sweep Inverse Filter SINV')
% ylabel('Amplitude (dB)')
% xlabel('Frequency (Hz)')

end

