clc; clear all; close all;

% Loads Data
[Kinetic1, Potential1, Energy1] = textread('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 2\Practice\ValentineCorrection\energy1.txt');
[Kinetic2, Potential2, Energy2] = textread('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 2\Practice\ValentineCorrection\energy2.txt');
[Kinetic3, Potential3, Energy3] = textread('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 2\Practice\ValentineCorrection\energy3.txt');
[Kinetic4, Potential4, Energy4] = textread('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 2\Practice\ValentineCorrection\energy4.txt');
[Kinetic5, Potential5, Energy5] = textread('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 2\Practice\ValentineCorrection\energy5.txt');
[Vx, Vy, Vz] = textread('C:\Users\bjcam\Google Drive\Molecular Dynamics\Project 2\Practice\8OclockCorrection\InitialVelocities.txt');

vel = [Vx, Vy, Vz];
velmean = mean(mean(vel'));

TimeStep = 1:1:1000;

% Energy
EnergyTrials = [Energy1, Energy2, Energy3, Energy4, Energy5];
EnergyMean = mean(EnergyTrials');
EnergySTD = std(EnergyTrials');

% Kinetic
KineticTrials = [Kinetic1, Kinetic2, Kinetic3, Kinetic4, Kinetic5];
KineticMean = mean(KineticTrials');
KineticSTD = std(KineticTrials');

% Potential
PotentialTrials = [Potential1, Potential2, Potential3, Potential4, Potential5];
PotentialMean = mean(PotentialTrials');
PotentialSTD = std(PotentialTrials');


figure;
set(gca,'fontsize',30)
errorbar(TimeStep, EnergyMean, EnergySTD, 'k');
ylim([85 105])
xlim([0 1000])
xlabel('Time Step', 'fontsize', 17);
ylabel('Reduced Energy', 'fontsize', 17);
title('Reduced Energy for Each Time-Step', 'fontsize', 17);


figure;
set(gca,'fontsize',30)
errorbar(TimeStep, KineticMean, KineticSTD, 'k');
%ylim([85 105])
xlim([0 1000])
xlabel('Time Step', 'fontsize', 17);
ylabel('Reduced Kinetic', 'fontsize', 17);
title('Reduced Kinetic Energy for Each Time-Step', 'fontsize', 17);

figure;
set(gca,'fontsize',30)
errorbar(TimeStep, PotentialMean, PotentialSTD, 'k');
%ylim([85 105])
xlim([0 1000])
xlabel('Time Step', 'fontsize', 17);
ylabel('Reduced Potential', 'fontsize', 17);
title('Reduced Potential Energy for Each Time-Step', 'fontsize', 17);
