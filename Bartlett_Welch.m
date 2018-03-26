% Bartlett's and Welch's Methods

clc
close all;
clear all;
fs=25e2; % Sampling frequency
Ts=1/fs; % Samplig period
f0=0.2*fs; % The frequency of the first sinusoid
f1=0.22*fs; % The frequency of the secon

%% Periodogram of two sinusoids

N=512; % Number of samples in the signal
t=0:Ts:(N-1)*Ts; % The time vector


% Generating 50 sample signals and their periodograms
for k=1:50
    phi= -pi + 2*pi.*rand(1,1);
    x(k,:)=15*sin(2*pi*f0*t+phi)+10*sin(2*pi*f1*t+phi)+randn(1,N);
    [X(k,:),w]=periodogram(x(k,:),[],2^12);
end
X_ave=mean(X);

figure;
% Overlay of the periodograms
X_dB=20*log10(X');
subplot(2,1,1);plot(w/(2*pi),X_dB)
axis([0 1/2 -80 10+max(max(X_dB))])
xlabel('Normalized Frequency')
ylabel('Magnitude Spectrum (dB)');
title('Periodogram')

% Average of the periodograms
X_ave_dB=20*log10(X_ave);
subplot(2,1,2);plot(w/(2*pi),X_ave_dB)
axis([0 1/2 -80 10+max(X_ave_dB)])

xlabel('Normalized Frequency')
ylabel('Magnitude Spectrum (dB)');

%% Bartlettt's method on the signal

N=512; % Number of samples in the signal
t=0:Ts:(N-1)*Ts; % The time vector
K=8; % Number of segments
L=N/K;
wind=ones(1,L); % Rectangular Window

% Generating 50 sample signals and their periodograms
for k=1:50
    phi= -pi + 2*pi.*rand(1,1);
    x(k,:)=15*sin(2*pi*f0*t+phi)+10*sin(2*pi*f1*t+phi)+randn(1,N);
    [X(k,:),w]=pwelch(x(k,:),wind,0,2^12);
end
X_ave=mean(X);

figure;
% Overlay of the periodograms
X_dB=20*log10(X');
subplot(2,1,1);plot(w/(2*pi),X_dB)
axis([0 1/2 -80 10+max(max(X_dB))])
title('Bartlett Method')
xlabel('Normalized Frequency')
ylabel('Magnitude Spectrum (dB)');

% Average of the periodograms
X_ave_dB=20*log10(X_ave);
subplot(2,1,2);plot(w/(2*pi),X_ave_dB)
axis([0 1/2 -80 10+max(X_ave_dB)])

xlabel('Normalized Frequency')
ylabel('Magnitude Spectrum (dB)');


%% Welch Method with 50% Overlap

N=512; % Number of samples in the signal
t=0:Ts:(N-1)*Ts; % The time vector
K=3; % Number of segments
L=2*N/(K+1);
wind=hamming(L); % Rectangular Window

% Generating 50 sample signals and their periodograms
for k=1:50
    phi= -pi + 2*pi.*rand(1,1);
    x(k,:)=15*sin(2*pi*f0*t+phi)+10*sin(2*pi*f1*t+phi)+randn(1,N);
    [X(k,:),w]=pwelch(x(k,:),wind,L/2,2^12);
end
X_ave=mean(X);

figure;
% Overlay of the periodograms
X_dB=20*log10(X');
subplot(2,1,1);plot(w/(2*pi),X_dB)
axis([0 1/2 -80 10+max(max(X_dB))])
title('Welch Method')
xlabel('Normalized Frequency')
ylabel('Magnitude Spectrum (dB)');

% Average of the periodograms
X_ave_dB=20*log10(X_ave);
subplot(2,1,2);plot(w/(2*pi),X_ave_dB)
axis([0 1/2 -80 10+max(X_ave_dB)])

xlabel('Normalized Frequency')
ylabel('Magnitude Spectrum (dB)');


