% This script 
% - picks a mat-file containing  Filtercoefficients in a HD-Filter object
% - picks a wav-file 
% - shows the frequeny- and  impulse_response of the filter
% - filters the wav file  using the filter-coefficients from the mat-file 
% - writes the filtered signal to wav_file
%...

clear all; 
close all;

% create inverse sweep
prompt= 'input speaker name';
ls_name=input(prompt)
prompt= 'length of the original sweep? (seconds)';
T = input(prompt); %sweep duration in secs (2^N)/fs;
prompt= 'start frequency of the original sweep? (Hz)';
f1=input(prompt); %starting frequency in Hz
prompt= 'end frequency of the original sweep? (Hz)';
f2=input(prompt); %ending frequency in Hz
fs=48000; %Samplingrate in Hz
Fs=48000; %Samplingrate in Hz

% loads filter and audio
[filename] = uigetfile('*.mat','Pick a filter-file');
filter_coefficients = load (filename);

[filename] = uigetfile('*.wav','Pick a wav-file');
[audio_in, fs] = audioread (filename);

[filename] = uigetfile('*.wav','Pick the test sweep (.wav)');
[test_sweep, fs] = audioread (filename);

% % generate Test-Sine 
%    dt = 1/fs;                   % seconds per sample
%    StopTime = 3;             % seconds
%    t = (0:dt:StopTime-dt)';     % seconds
%    %%Sine wave:
%    Fc = 60;                     % hertz
%    audio_in = 0.5*sin(2*pi*Fc*t);

%%%% fade in audio!
ft_samples=15;
fout_samples=15;
fade_vector=linspace(0,1,ft_samples);
fout_vector=linspace(1,0,fout_samples);
begin_fout=length(audio_in)-length(fout_vector)+1;

%%%% filters audio!
audio_in_fade=audio_in;
audio_in_fade(1:ft_samples)=audio_in(1:ft_samples).*fade_vector';
audio_in_fade(begin_fout:length(audio_in_fade))=audio_in_fade(begin_fout:length(audio_in_fade)).*fout_vector';
filtered_audio = filter(filter_coefficients.Hd, audio_in_fade);

%%%% PLOT FILTERED AUDIO

%%%% Ignore
% axis for Samples
% a=axes('Position',[.1 .1 .8 1e-12]);
%set(a,'Units','normalized');
%set(a,'Color','none');
% axis for Frequency
% b=axes('Position',[.1 .2 .8 .7]);
% set(b,'Units','normalized');
% log2end=f1+((f2-f1)/T)*(length(filtered_audio)/fs)
% % set limits and labels
% set(a,'xlim',[f1 (log2end)]);
% set(b,'xlim',[1 T*fs]);
% xlabel(a,'Frequency (not exact)')
% xlabel(b,'Samples')

figure(1)
figure('name','Compare filtered and unfiltered signal (time domain)')
subplot(3,1,1)
plot(filtered_audio,'r');
title('filtered and faded audio')
ylim([-max(filtered_audio)-0.01,max(filtered_audio)+0.01])
xlabel('samples')
hold on

subplot(3,1,2)
plot(audio_in,'b');
title('recorded signal')
ylim([-max(filtered_audio),max(filtered_audio)])
xlabel('samples')


subplot(3,1,3)
plot(test_sweep,'g');
title('original test sweep')
ylim([-max(test_sweep)-0.3,max(test_sweep)+0.3])
xlabel('samples')
hold off

%%%% Show filter details
freqz(filter_coefficients.Hd);
pause(2);
%%% show impulse-response
impz(filter_coefficients.Hd);
pause(2)

%%%%% Compensation Filterdelay for FIR Filter
%defines len_filter according to Numerator of the chosen filter 
% len_filter = length(filter_coefficients.Hd.Numerator);
% new_zeroes = zeros(1,floor(len_filter/2));
% audio_in_extend = horzcat(audio_in',new_zeroes);
% filtered_audio = filter(filter_coefficients.Hd, audio_in_extend);

%%%% Compensation Filterdelay (FIR) Filtered Signal
% concat_signal = horzcat(new_zeroes,audio_in');
%plot(concat_signal,'b');
  
% Audio wird in Spektralbereich gewandelt;
audio_in_fft = fft(audio_in);
Q2 = abs(audio_in_fft/length(audio_in));
Q1 = Q2(1:length(audio_in)/2+1);
Q1(2:end-1) = 2*Q1(2:end-1);

filtered_audio_fft = fft(filtered_audio);
P2 = abs(filtered_audio_fft/length(filtered_audio));
P1 = P2(1:length(filtered_audio)/2+1);
P1(2:end-1) = 2*P1(2:end-1);

testsweep_fft = fft(test_sweep);
R2 = abs(testsweep_fft/length(test_sweep));
R1 = R2(1:length(test_sweep)/2+1);
R1(2:end-1) = 2*R1(2:end-1);

f = Fs*(0:(length(audio_in)/2))/length(audio_in);
g = Fs*(0:(length(filtered_audio)/2))/length(filtered_audio);
h = Fs*(0:(length(test_sweep)/2))/length(test_sweep);

figure(3)
subplot(3,1,1)

semilogx(f,Q1, 'b')
title('Frequency Spectrum of the unfiltered exp. Sine Sweep')
ylabel('Amplitude (dB)')
xlabel('Frequency (Hz)')

figure(3)
subplot(3,1,2)
semilogx(g,P1, 'r')
title('Amplitude Spectrum of the filtered exp. Sine Sweep SIG')
ylabel('Amplitude (dB)')
xlabel('Frequency (Hz)')

figure(3)
subplot(3,1,3)
semilogx(h,R1, 'g')
title('Amplitude Spectrum of the original test_sweep')
ylabel('Amplitude (dB)')
xlabel('Frequency (Hz)')

% Create the swept sine tone
w1 = 2*pi*f1; %Omega1
w2 = 2*pi*f2; %Omega2
K = T*w1/log(w2/w1); %growingparameter for creating the sweep
L = T/log(w2/w1); %lengthparameter for creating the sweep
t = linspace(0,T-1/fs,fs*T); %time vector
s = sin(K*(exp(t/L) - 1)); %exponential sweep timesignal
% LT=length(s); %length of sweep vector

% find max value for filtered signal to normalize invsweep in time dimension.
maximum = max(filtered_audio);

invsweep=s(length(s):-1:1).*exp(-t./L);%generate backwards exponential sine sweep
% sinv=invsweep/max(abs(invsweep)); %normalize timesignal to maximum value of time signal
sinv=invsweep*maximum; %normalize timesignal to maximum value of time signal

sigL=filtered_audio';
zero= zeros(1,length(filtered_audio)-length(invsweep)); %% zeropadding 
sinvL=horzcat(sinv,zero);
NFFT=length(sigL); %set FFT window length to sweep length

SIG=fft(sigL);
SINV=fft(sinvL);

figure(4) 

% PLOT SIG AND SINV (spectral)
f = fs/2*linspace(0,1,NFFT/2); %create plotting vector
subplot(3,1,1), loglog(f,abs(SIG(1:NFFT/2)/max(abs(SIG)))) %Plot single-sidedamplitude spectrum of the sweep
title('Amplitude Spectrum of the filtered exp. Sine Sweep SIG')
ylabel('Amplitude (dB)')
xlabel('Frequency (Hz)')
subplot(3,1,2), loglog(f,abs(SINV(1:NFFT/2))/max(abs(SINV))) %Plot single-sidedamplitude spectrum of inverse filter
title('Amplitude Spectrum of the Sine Sweep Inverse Filter SINV')
ylabel('Amplitude (dB)')
xlabel('Frequency (Hz)')

%DECONVOLVE Sweep with corrected inverse filter 
IMP=SIG.*SINV;
imp=ifft(IMP);
imp=imp.*max(imp);

subplot(3,1,3)
plot(20*log10(abs(imp/max(imp))))
title('recovered impulse response')
ylabel('Amplitude (dB)')
xlabel('Time (samples)')
axis([NFFT*0.7 NFFT*1.1 -110 0]);

figure(5)
subplot(2,1,1)
% semilogx(f,(20*log10(abs(IMP(1:NFFT/2))/max(abs(IMP)))))
% loglog(f,((abs(IMP(1:NFFT/2))/max(abs(IMP)))))
norm1k=abs([IMP(round(1000/(fs/NFFT)))]);
% semilogx(f,20*log10(abs(IMP(1:NFFT/2)/max(abs(IMP))))) %normalize to
% global max
normlog=20*log10(abs(IMP(1:NFFT/2)/norm1k));
semilogx(f,normlog); %normalize to 1kHz

grid on
title('amplitude-frequency response (normalized to 1kHz)')
ylabel('Amplitude (dB)')
xlabel('f (Hz)')
xlim([40 22000])
ylim([-12,8])

%mittelwertsbildung der intervalle evtl ueberlappung

% test sweep kill first 460 samples
% testsweep delete
% optisch prüfen 460 samples
% abschneiden 460 samples
% +2 peridoen abschneiden beibeiden signalen
% 500
% neues testsignal generieren mit Anfang = 1. periodendauer
% zusätzliche verschiebung
% 924 = T = 51,948 Hz at 48k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=test_sweep;
s=s';
y1=filtered_audio;
diff=(length(s)-length(y1));
zeroblock = zeros(1,diff);
y=horzcat(y1',zeroblock);

LT=length(s);

nfft=length(s);

xPlot = (linspace(0,fs/2,((nfft/2)+1))); %create plotting vector

window=rectwin(nfft);
noverlap=0;
x=s;

[Txy,F] = tfestimate(x,y,window,noverlap,nfft,fs);

figure(16)
plot(xPlot,angle(Txy));
title('angle wraped in rads 1rad=3,14=180deg');
hold off

tol = pi;

figure(17);
plot(xPlot,(unwrap(angle(Txy),tol)));
title('angle unwraped 1rad=3,14=180deg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot impulse response in time domain
figure(7)
[impmax, Samplemax] = max(imp);
Timpmax = (Samplemax/fs);
tSampling = 1/fs;
n=length(imp);
t=0:tSampling:0+(n-1)*tSampling;
impnorm=20*log10(abs(imp./impmax)); %normalize imp
plot((t-Timpmax),impnorm);
xticks(-0.01:0.01:0.05)
xlim([-0.03,0.05])
% ylim([-1*max(impnorm)/1.25,max(impnorm)/1.25])
ylabel('normalized amplitude (linear)')
xlabel('t (seconds)')
title('aligned impulse response (time domain)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
