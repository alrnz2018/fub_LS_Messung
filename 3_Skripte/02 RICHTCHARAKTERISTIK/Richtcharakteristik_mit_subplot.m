close all;
clear;

addpath('./Messdaten');
%% Import the data
% Erstelle Matrix f�r alle Messdaten mit 46 Winkelpositionen f�r zwei
% Mikrofone; 148000 Zeilen damit alle Dateien prolemlos ausgelesen werden
import_data = zeros(148000,180);

b = 1;
for (a = 0:2:90)
    a1 = int2str(a);
    Mic1 = strcat('./Messdaten/',a1,' Grad/Messung_ch1'," ",a1,' Grad.wav');
    Mic2 = strcat('./Messdaten/',a1,' Grad/Messung_ch2'," ",a1,' Grad.wav');
    [read1, FS] = audioread(Mic1);
    [read2, FS] = audioread(Mic2);
    mic1tospalte=b;
    %mic1tospiegel=181-b;
    mic2tospalte=93-b;
    import_data(1:1:length(read1),mic1tospalte) = read1;
    %import_data(1:1:length(read1),mic1tospiegel) = read1;
    import_data(1:1:length(read2),mic2tospalte) = read2;
    b = b + 1;
end
%Verzichte auf Mittelung des 90 Grad Wertes, weil Spalte 46 und 47 selber
%messpunkt

%Definiere kompletten Kreis
Spiegel_mat = zeros (148000,180);

Spiegel_mat(:,1) = import_data(:,1);
Spiegel_mat(:,91) = import_data(:,92);

for (a = 2:1:90)
    Spiegel_mat(:,a) = import_data(:,a);
    b = 182 - a;
    Spiegel_mat(:,b) = import_data(:,a);
end


%cd Messdaten
%dir *.wav

%% FFT

Fourier_mat = fft(Spiegel_mat);

%km184 = cell2mat(raw);

%km184_f = km184(:,1);
%km184_P = km184;
%km184_P (:,1) = [];
%%

%km184_p = 10.^(km184_P/20);   % von dB in Druck

% Normieren auf 0°
%for m = 1:length(km184_p(:,1))
%km184_p(m,:) = km184_p(m,:)./km184_p(m,1);
%end

% Winkel-Vektor berechnen (0-360°)
rad = (0:2:358) .* (pi / 180);

% Frequenzvektor in Oktaven von 125Hz bis 16kHz
freq = [125, 250, 489, 980, 1960, 3945, 7910, 1586];
% freq = [125, 250, 500, 1000, 2000, 4000, 7910, 16000];


%%

% Datenreihe von 0° - 360° für die Oktaven, normiert auf 0 Grad
%km184_oktav = oktavread( km184_p, km184_f, freq);


% Plot der Polardiagramme
%g5 = figure;
%%mmpolar(rad,km184_oktav(1,:),rad,km184_oktav(2,:),rad,km184_oktav(3,:),...
%        rad,km184_oktav(4,:),rad,km184_oktav(5,:),rad,km184_oktav(6,:),...
%        rad,km184_oktav(7,:),rad,km184_oktav(8,:), 'Style', 'compass', 'RLimit', ...
%        [-25 5], 'RTickValue', [5 0 -5 -10 -15 -20 -25], 'TTickDelta', 45, freq);
%title('Polardiagramm KM184');
%legend('125Hz','250Hz','500Hz','1000Hz','2000Hz','4000Hz','8000Hz','16000Hz','Location','BestOutside');
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





