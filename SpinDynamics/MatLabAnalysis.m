clear all; clc; close all;

L10Trail1 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10Trail1.txt');
L10Trail2 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10Trail2.txt');
L10Trail3 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10Trail3.txt');
L10Trail4 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10Trail4.txt');
L10Trail5 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10Trail5.txt');

L24Trail1 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL24Trail1.txt');
L24Trail2 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL24Trail2.txt');
L24Trail3 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL24Trail3.txt');
L24Trail4 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL24Trail4.txt');
L24Trail5 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL24Trail5.txt');

L74Trail1 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL74Trail1.txt');
L74Trail2 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL74Trail2.txt');
L74Trail3 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL74Trail3.txt');
L74Trail4 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL74Trail4.txt');
L74Trail5 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL74Trail5.txt');


L10T01Trail1 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10T01Trail1.txt');
L10T01Trail2 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10T01Trail2.txt');
L10T01Trail3 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10T01Trail3.txt');
L10T01Trail4 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10T01Trail4.txt');
L10T01Trail5 = load('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 5\Data\CSVOutputs\CorrL10T01Trail5.txt');

time = 1:1:500;
fftTime = linspace(0,2*pi, 500);

time2 = 1:1:2500;
fftTime2 = linspace(0,2*pi, 2500);



%% Correlation L10
CorrL10 = [L10Trail1, L10Trail2, L10Trail3, L10Trail4, L10Trail5];
MeanL10 = mean(CorrL10');
STDL10 = std(CorrL10');

FFTl10 = fft(MeanL10);

figure;
subplot(2,1,1);
errorbar(time, MeanL10, STDL10);
title('Correlation Function for L = 10, TimeStep = 0.1');

xlabel('Time Separation (J^-^1)');
ylabel('Correlation');
set(gca,'fontsize',17)

subplot(2,1,2);
plot(fftTime, abs(FFTl10));
title('Dynamic Structure Factor for L = 10, TimeStep = 0.1');
xlabel('\omega J^-^1)')
ylabel('S(q,\omega)')
xlim([0 2*pi])
set(gca,'fontsize',17)



%% Correlation L24
CorrL24 = [L24Trail1, L24Trail2, L24Trail3, L24Trail4, L24Trail5];
MeanL24 = mean(CorrL24');
STDL24 = std(CorrL24');
FFTl24 = fft(MeanL24);

figure;
subplot(2,1,1);
errorbar(time, MeanL24, STDL24);
title('Correlation Function for L = 24, TimeStep = 0.1');
xlabel('Time Separation (J^-^1)');
ylabel('Correlation');
set(gca,'fontsize',17)


subplot(2,1,2);
plot(fftTime, abs(FFTl24));
title('Dynamic Structure Factor for L = 24, TimeStep = 0.1');
xlabel('\omega J^-^1')
ylabel('S(q,\omega)')
xlim([0 2*pi])
set(gca,'fontsize',17)


%% Correlation L74
CorrL74 = [L74Trail1, L74Trail2, L74Trail3, L74Trail4, L74Trail5];
MeanL74 = mean(CorrL74');
STDL74 = std(CorrL74');
FFTl74 = fft(MeanL74);


figure;
subplot(2,1,1);
errorbar(time, MeanL74, STDL74);
title('Correlation Function for L = 74, TimeStep = 0.1');
xlabel('Time Separation (J^-^1)');
ylabel('Correlation');
set(gca,'fontsize',17)

subplot(2,1,2);
plot(fftTime, abs(FFTl74));
title('Dynamic Structure Factor for L = 74, TimeStep = 0.1');
xlabel('\omega J^-^1')
ylabel('S(q,\omega)')
xlim([0 2*pi])
set(gca,'fontsize',17)

%% Correlation L10 TimeStep 0.01
CorrL10T01 = [L10T01Trail1, L10T01Trail2, L10T01Trail3, L10T01Trail4, L10T01Trail5];
MeanL10T01 = mean(CorrL10T01');
STDLL10T01 = std(CorrL10T01');
FFTlL10T01 = fft(MeanL10T01);

figure;
subplot(2,1,1);
errorbar(time2, MeanL10T01, STDLL10T01);
title('Correlation Function for L = 10, TimeStep = 0.01');
xlabel('Time Separation (J^-^1)');
ylabel('Correlation');
set(gca,'fontsize',17)

subplot(2,1,2);
plot(fftTime2, abs(FFTlL10T01));
title('Dynamic Structure Factor for L = 10, TimeStep = 0.01');
xlabel('\omega J^-^1')
ylabel('S(q,\omega)')
xlim([0 2*pi])
set(gca,'fontsize',17)