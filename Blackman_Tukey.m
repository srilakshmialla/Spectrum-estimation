% EE-658 Advanced Digital Signal Processing
% Blackman-Tukey Method

clc
close all;
clear all;

fs=25e2; % Sampling frequency
Ts=1/fs; % Samplig period
f0=0.2*fs; % The frequency of the first sinusoid
f1=0.22*fs; % The frequency of the secon

%% Periodogram of two sinusoids

N=512; % Number of samples in the signal
M=400;
t=0:Ts:(N-1)*Ts; % The time vector
wind=hamming(2*M-1);

% Generating 50 sample signals and their periodograms
for k=1:50
    phi= -pi + 2*pi.*rand(1,1);
    x(k,:)=15*sin(2*pi*f0*t+phi)+10*sin(2*pi*f1*t+phi)+randn(1,N);
    [X(k,:),w]=per_smooth(x(k,:),wind,M);
end
X_ave=mean(X);

figure;
% Overlay of the periodograms
X_dB=20*log10(X');
subplot(2,1,1);plot(w,X_dB)
axis([0 pi -80 10+max(max(X_dB))])
xlabel('Normalized Frequency')
ylabel('Magnitude Spectrum (dB)');
title('B-T Spectrum')

% Average of the periodograms
X_ave_dB=20*log10(X_ave);
subplot(2,1,2);plot(w,X_ave_dB)
axis([0 pi -80 10+max(X_ave_dB)])

xlabel('Normalized Frequency')
ylabel('Magnitude Spectrum (dB)');