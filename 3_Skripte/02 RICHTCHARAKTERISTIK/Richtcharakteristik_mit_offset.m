%% Richtcharakteristikskript
% Author: Tobias Festag, Adrian Lorenz, August 2019
% Plottet aus 46 Messdatensätzen á 2 Signalen (um 180 Grad versetzt) ein Polardiagramm unterteilt in Oktavbänder.
% Das Messsignal ist ein logarithmischer Sweep, eine entsprechende Kompensation wird später vorgenommen.


close all;
clear all;

%% Lade Störsignalfilter (FIR)
[filename] = uigetfile('*.mat','Pick a filter-file');
filter_coefficients = load (filename);
filter_coefficients = load ('E:\20170109_LS_Messung_AuswertungSS19\02 RICHTCHARAKTERISTIK\Filter\FIR_4499_Order_Window_Hann_-6dB_60Hz_FS48Khz.mat');

%% Lade Messdatenverzeichnis
prompt= 'Name your Speaker: ';
Speaker = input(prompt,'s'); 
[path_speaker] = uigetdir('./','Select your Speakerdata');
addpath(path_speaker);
%% Import the data
% Erstelle Matrix für alle Messdaten mit 46 Winkelpositionen für zwei
% Mikrofone; 148000 Zeilen damit alle Dateien prolemlos ausgelesen werden

import_data = zeros(148000,180);

b = 1;
for a = 0:2:90
    a1 = int2str(a);
    Mic1 = strcat(path_speaker,'/',a1,' Grad/Messung_ch1'," ",a1,' Grad.wav');
    Mic2 = strcat(path_speaker,'/',a1,' Grad/Messung_ch2'," ",a1,' Grad.wav');
    [read1, FS] = audioread(Mic1);
    [read2, FS] = audioread(Mic2);
        
    %% filtere Störsignal raus und fade das Nutzsignal
    % Kompensiere Lautfzeit des FIR
    % noch fehlerhaft
    filtered1 = fade_and_filter(read1,filter_coefficients);
    filtered2 = fade_and_filter(read2,filter_coefficients);

    mic1tospalte=b;
    mic2tospalte=93-b;

    import_data(1:1:length(filtered1),mic1tospalte) = filtered1;
    import_data(1:1:length(filtered2),mic2tospalte) = filtered2;

    b = b + 1;
end
%% Aufnahmeseitigen Offsetfaktor zwischen Mic1 und Mic2? Nach Importschleife stehen die Importvariablen beim gleichen Wert
% offset = (sum(filtered2)-sum(filtered1))/length(filtered1);
% diff_vektor = filtered2 - filtered1;
% offset = sum(diff_vektor)/length(diff_vektor);

offset = 1;

%Definiere kompletten Kreis
Spiegel_mat = zeros (length(import_data),180);

%offset = vertcat(offset,zeros((size(Spiegel_mat,1)-length(offset)),1)) % Prolongue offset vector length of array

Spiegel_mat(:,1) = import_data(:,1) * offset;
Spiegel_mat(:,91) = import_data(:,92);

for a = 2:1:90
	% Multipliziere Offsetfaktor mit Messdaten von 0°-90° und 270° - 360° zum Ausgleich etwaiger Aussteuerungsungleichheiten
	if a<=46 
    	Spiegel_mat(:,a) = import_data(:,a) * offset;
   		b = 182 - a;
    	Spiegel_mat(:,b) = import_data(:,a) * offset;
    else
    	Spiegel_mat(:,a) = import_data(:,a);
   		b = 182 - a;
    	Spiegel_mat(:,b) = import_data(:,a);
    end
end

%% FFT unserer Messignalmatrix

Fourier_mat = fft(Spiegel_mat);
Abs_mat = abs(Fourier_mat);


%%
% generiere inversen sweep

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

name_plot = "Polardiagramm " + Speaker;
 title(name_plot);
 legend('125Hz','250Hz','500Hz','1000Hz','2000Hz','4000Hz','8000Hz','16000Hz','Location','BestOutside');
 % thetalim([-90 90]);
 