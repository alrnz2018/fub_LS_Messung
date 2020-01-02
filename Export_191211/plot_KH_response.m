close all

fs=48000;
NFFT=145674;
f = fs/2*linspace(0,1,NFFT/2);
f = f(1:72812);

figure(1)

KH120_KC1=KH120_KC1(1:72812);
KH120_KC2=KH120_KC2(1:72812);
KH120_MG1=KH120_MG1(1:72812);

semilogx(f,KH120_KC1);
hold on
semilogx(f,KH120_KC2);
hold on
semilogx(f,KH120_MG1);
hold on
semilogx(f,KH120_MG2); 

grid on
title('amplitude-frequency response (normalized to 1kHz)')
ylabel('Amplitude (dB)')
xlabel('f (Hz)')
xlim([40 22000])
ylim([-12,8])

legend('KH120-KC1','KH120-KC2','KH120-MG1','KH120-MG2');

figure(2)

order = 9;
framelen = 29;

smoo_KH120_KC1=sgolayfilt(KH120_KC1,order,framelen);
smoo_KH120_KC2=sgolayfilt(KH120_KC2,order,framelen);
smoo_KH120_MG1=sgolayfilt(KH120_MG1,order,framelen);
smoo_KH120_MG2=sgolayfilt(KH120_MG2,order,framelen);

semilogx(f,smoo_KH120_KC1);
hold on
semilogx(f,smoo_KH120_KC2);
hold on
semilogx(f,smoo_KH120_MG1);
hold on
semilogx(f,smoo_KH120_MG2);
hold on

grid on
title('sgolayfilt order 3 framelength 11 amplitude-frequency response (normalized to 1kHz)')
ylabel('Amplitude (dB)')
xlabel('f (Hz)')
xlim([40 22000])

ylim([-40,6])
legend('KH120-KC1','KH120-KC2','KH120-MG1','KH120-MG2');

