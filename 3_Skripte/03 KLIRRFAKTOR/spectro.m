clear all;
close all;


res = 1020; %Aufloesung des Spektrogrammes uebliche Werte zwischen 400 (grob) und MAX (!!!) 1023 (fein)

%%Matlab Skript zur Erstellung eines Spektrogrammes einer Audio Datei
%
%Author: Tobias Festag
%
%Built: 14.08.2016
%Nicht Benutzerfreundliche Oberflaeche! Skript muss im gleichen ordner
%%ausgefuehrt werden in dem die zu analysierende datei mit Namen '01.wav', liegt. Aufloesung muss im Quellcode eingestellt werden 

disp('Willkommen!')
[filename] = uigetfile('*.wav','Pick a sound file (WAV)');
[audio, fs] = audioread(filename);
[X,F,T,P] = spectrogram(audio, 1024, res);
surf(T,F,10*log10(abs(P)),'EdgeColor','none');