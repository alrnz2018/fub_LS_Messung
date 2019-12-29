clear all;
close all;

% Abfrage wie viele Files eingelesen werden sollen

prompt = 'Wie viele Files sollen eingelesen werden: ';
b= input(prompt);

%Lese Audiodateien ein, füge die Pfade zum Workspace hinzu und plotte die Spektrogramme
for a=1:1:b

	prompt = 'Benennen Sie das Signal eindeutig: ';
	sig_name = input(prompt, 's');
	[file,path] = uigetfile('*.wav','Pick a signal file');
	addpath(path);

	[audio, fs] = audioread(file);

	figure;
	s = spectrogram(audio);
	spectrogram(audio,'yaxis')
	title(sig_name);

end
