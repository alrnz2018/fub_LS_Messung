close all;
clear all;

%% Lade Störsignalfilter (FIR)
[filename] = uigetfile('*.mat','Pick a filter-file');
%filter_coefficients = load (filename);
filter_coefficients = load ('E:\20170109_LS_Messung_AuswertungSS19\02 RICHTCHARAKTERISTIK\Filter\FIR_4499_Order_Window_Hann_-6dB_60Hz_FS48Khz.mat');

%%
addpath('./Messdaten_KH120');
%% Import the data
% Erstelle Matrix für alle Messdaten mit 46 Winkelpositionen für zwei
% Mikrofone; 148000 Zeilen damit alle Dateien prolemlos ausgelesen werden
import_data = zeros(148000,180);

b = 1;
for (a = 0:2:90)
    a1 = int2str(a);
    Mic1 = strcat('./Messdaten_KH120/',a1,' Grad/Messung_ch1'," ",a1,' Grad.wav');
    Mic2 = strcat('./Messdaten_KH120/',a1,' Grad/Messung_ch2'," ",a1,' Grad.wav');
    [read1, FS] = audioread(Mic1);
    [read2, FS] = audioread(Mic2);
    
    %% FILTER AND FADE SPÄTER IN SEPARATE FUNKTION
    
    
    
    %% filtere Störsignal raus und fade das Nutzsignal
    % Kompensiere Lautfzeit des FIR
    % noch fehlerhaft
    filtered1 = fade_and_filter(read1,filter_coefficients);
    filtered2 = fade_and_filter(read2,filter_coefficients);
    
    %%
%   entweder    
    mic1tospalte=b;
    mic2tospalte=93-b;
    %     mic1tospiegel=181-b;  ---> eigentlich nutzlos!!!
%   oder
%     mic1tospiegel=b;
%     mic2tospiegel=93-b;
    % entweder
    import_data(1:1:length(filtered1),mic1tospalte) = filtered1;
    import_data(1:1:length(filtered2),mic2tospalte) = filtered2;
    % oder
%     import_data(1:1:length(read1),mic1tospiegel) = read1;    
%     import_data(1:1:length(read2),mic2tospiegel) = read2;
    b = b + 1;
end
%Aufnahmeseitigen Offsetfaktor Mic1 zu Mic2?
offset = filtered2 - filtered1;

%Definiere kompletten Kreis
Spiegel_mat = zeros (148000,180);

Spiegel_mat(:,1) = import_data(:,1);
Spiegel_mat(:,91) = import_data(:,92);

for (a = 2:1:90)
    Spiegel_mat(:,a) = import_data(:,a);
    b = 182 - a;
    Spiegel_mat(:,b) = import_data(:,a);
end

%% FFT unserer Messignalmatrix

Fourier_mat = fft(Spiegel_mat);
Abs_mat = abs(Fourier_mat);


%%
%generiere inversen sweep

prompt= 'length of the original sweep? (seconds)';
T = input(prompt); %sweep duration in secs (2^N)/fs;
prompt= 'start frequency of the original sweep? (Hz)';
f1=input(prompt); %starting frequency in Hz
prompt= 'end frequency of the original sweep? (Hz)';
f2=input(prompt); %ending frequency in Hz
max_sig = max(max(Abs_mat)); 
sinv=gen_invsweep(T,f1,f2,FS,max_sig); % generiere Sweep und normalisiere ihn auf Maximalwert der Matrix


zero= zeros(1,length(Abs_mat)-length(sinv)); 
sinv=horzcat(sinv,zero); % adjust length of invsweep to length of testsignal
abs_SINV=abs(fft(sinv)).*max_sig; %normalized sinv fft

UnPink_mat=Abs_mat.*abs_SINV';

%%
%cd Messdaten
%dir *.wav

%%
% Oktavbandvektor bestimmen
freq = [125, 250, 500, 1000, 2000, 4000, 7910, 16000];
fft_resolution=FS/length(UnPink_mat);

freq_lower = floor(freq./1.41);
freq_upper = floor(freq./0.7);

for i = 1:length(freq)
    
N_lower(i) = floor(freq_lower(i)/fft_resolution);
N_upper(i) = floor(freq_upper(i)/fft_resolution);
N_width(i) = N_upper(i)-N_lower(i);

for j = 1:1:180
    Temp_value = UnPink_mat(N_lower(i):N_upper(i),j);
    Abs_Octa(i,j) = sum(Temp_value)/N_width(i);
end
end


%% Plotte

%km184_p = 10.^(km184_P/20);   % von dB in Druck

% Normieren auf 0Â°
%for m = 1:length(km184_p(:,1))
%km184_p(m,:) = km184_p(m,:)./km184_p(m,1);
%end

% Winkel-Vektor berechnen (0-360Â°)
rad = (0:2:358) .* (pi / 180);

% Frequenzvektor in Oktaven von 125Hz bis 16kHz
% freq = [125, 250, 489, 980, 1960, 3945, 7910, 1586];
% freq = [125, 250, 500, 1000, 2000, 4000, 7910, 16000];


%%

% Datenreihe von 0Â° - 360Â° fÃ¼r die Oktaven, normiert auf 0 Grad
%km184_oktav = oktavread( km184_p, km184_f, freq);


%% Plot der Polardiagramme
g = figure;
polarplot(rad,Abs_Octa(1,:),rad,Abs_Octa(2,:),rad,Abs_Octa(3,:),...
         rad,Abs_Octa(4,:),rad,Abs_Octa(5,:),rad,Abs_Octa(6,:),...
         rad,Abs_Octa(7,:),rad,Abs_Octa(8,:));
% Troubleshooting Backup
% mmpolar(rad,Abs_Octa(1,:),rad,Abs_Octa(2,:),rad,Abs_Octa(3,:),...
%          rad,Abs_Octa(4,:),rad,Abs_Octa(5,:),rad,Abs_Octa(6,:),...
%          rad,Abs_Octa(7,:),rad,Abs_Octa(8,:), 'Style', 'compass', 'RLimit', ...
%          [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);

 title('Polardiagramm KH120');
 legend('125Hz','250Hz','500Hz','1000Hz','2000Hz','4000Hz','8000Hz','16000Hz','Location','BestOutside');
 % thetalim([-90 90]);
 
% Speichern als PDF
%save_plot(g5, 'nr4a_polar');

%g6 = figure;
%subplot(2,2,1)
%mmpolar(rad,km184_oktav(1,:), 'Style', 'compass', 'RLimit', [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);
%title('Polardiagramm KM184 bei 125 Hz');

%subplot(2,2,2)
%mmpolar(rad,km184_oktav(2,:), 'Style', 'compass', 'RLimit', [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);
%title('Polardiagramm KM184 bei 250 Hz');

%subplot(2,2,3)
%mmpolar(rad,km184_oktav(3,:), 'Style', 'compass', 'RLimit', [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);
%title('Polardiagramm KM184 bei 500 Hz');

%subplot(2,2,4)
%mmpolar(rad,km184_oktav(4,:), 'Style', 'compass', 'RLimit', [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);
%title('Polardiagramm KM184 bei 1000 Hz');

% Speichern als PDF
%save_plot(g6, 'nr4b_polar');

%g7 = figure;
%subplot(2,2,1)
%mmpolar(rad,km184_oktav(5,:), 'Style', 'compass', 'RLimit', [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);
%title('Polardiagramm KM184 bei 2000 Hz');

%subplot(2,2,2)
%mmpolar(rad,km184_oktav(6,:), 'Style', 'compass', 'RLimit', [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);
%title('Polardiagramm KM184 bei 4000 Hz');

%subplot(2,2,3)
%mmpolar(rad,km184_oktav(7,:), 'Style', 'compass', 'RLimit', [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);
%title('Polardiagramm KM184 bei 8000 Hz');

%subplot(2,2,4)
%mmpolar(rad,km184_oktav(8,:), 'Style', 'compass', 'RLimit', [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);
%title('Polardiagramm KM184 bei 16000 Hz');

% Speichern als PDF
%save_plot(g7, 'nr4c_polar');





