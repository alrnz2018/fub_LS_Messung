function [filtered_audio] = fade_and_filter(audio_in,filter_coefficients)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% FIR Fitler Designen mit Stopband unter 50 Hz
% FIR Filter als Variable fest einbauen, damit mans nicht neu einladen muss bei jedem
% Durchlauf

% fades audio!
ft_samples=15; %fadein sample
fout_samples=15;
fade_vector=linspace(0,1,ft_samples);
fout_vector=linspace(1,0,fout_samples);
begin_fout=length(audio_in)-length(fout_vector)+1;
audio_in_fade=audio_in;
audio_in_fade(1:ft_samples)=audio_in(1:ft_samples).*fade_vector';
audio_in_fade(begin_fout:length(audio_in_fade))=audio_in_fade(begin_fout:length(audio_in_fade)).*fout_vector';

% Compenstion Filterdelay for FIR Filter
% defines len_filter according to Numerator of the chosen filter 

len_filter = length(filter_coefficients.Hd.Numerator);
new_zeroes = zeros(1,floor(len_filter/2));
audio_in_extend = horzcat(audio_in_fade',new_zeroes);

% filters audio!
filtered_audio = filter(filter_coefficients.Hd, audio_in_extend);

% Compensation Filterdelay (FIR) Filtered Signal
%concat_signal = audio_in'((new_zeroes+1:);
len=length(filtered_audio);
concat_signal = filtered_audio(length(new_zeroes)+1:length(filtered_audio));
filtered_audio=concat_signal';

% Vektorbeschneidung soll am ende wieder bei 145540 Werten liegen! zeile 32
% macht aber noch nicht was es soll
% überprüfe Matrizenausrichtung der Vektoren für die Wertübergabe
end

