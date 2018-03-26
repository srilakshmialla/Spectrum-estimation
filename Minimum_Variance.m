% Minimum Variance

clc
close all;
clear all;

fs=25e2; % Sampling frequency
Ts=1/fs; % Samplig period
f0=0.2*fs; % The frequency of the first sinusoid
f1=0.22*fs; % The frequency of the secon


N=512; % Number of samples in the signal
p=512;
t=0:Ts:(N-1)*Ts; % The time vector

% Generating 50 sample signals and their periodograms
phi= -pi + 2*pi.*rand(1,1);
x=15*sin(2*pi*f0*t+phi)+10*sin(2*pi*f1*t+phi)+randn(1,N);
Nf=2^11; % Number of samples in the FFT
f=0:fs/Nf:(Nf-1)/Nf*fs/2;
X=MV(x,p,f);

figure;
% Overlay of the periodograms
X_dB=10*log10(X);
plot(f,X_dB)
axis([0 fs/2 min(min(X_dB))-10 10+max(max(X_dB))])
xlabel('Normalized Frequency')
ylabel('Magnitude Spectrum (dB)');
title('M-V Spectrum')
