% EE-658 Advanced Digital Signal Processing
% Periodogram Examples


%% Random phase sinusoid in the white noise
clc
close all;
clear;

fs=25e2; % Sampling frequency
Ts=1/fs; % Samplig period
f0=0.2*fs; % The frequency of the signal

N=64; % Number of samples in the signal
t=0:Ts:(N-1)*Ts; % The time vector


% Generating 50 sample signals and their periodograms
for k=1:50
    phi= -pi + 2*pi.*rand(1,1);
    x(k,:)=5*sin(2*pi*f0*t+phi)+randn(1,N);
    [X(k,:),w]=periodogram(x(k,:),[],2^12);
end
X_ave=mean(X);

% Overlay of the periodograms
X_dB=20*log10(X');
figure;plot(w/(2*pi),X_dB)
axis([0 1/2 -80 10+max(max(X_dB))])

% Average of the periodograms
X_ave_dB=20*log10(X_ave);
figure;plot(w/(2*pi),X_ave_dB)
axis([0 1/2 -80 10+max(X_ave_dB)])

%% Two random phase sinusoid in the white noise

clc
clear all;

fs=25e2; % Sampling frequency
Ts=1/fs; % Samplig period
f0=0.2*fs; % The frequency of the first sinusoid
f1=0.21*fs; % The frequency of the second sinusoid

N=256; % Number of samples in the signal
t=0:Ts:(N-1)*Ts; % The time vector


% Generating 50 sample signals and their periodograms
for k=1:50
    phi= -pi + 2*pi.*rand(1,1);
    x(k,:)=5*sin(2*pi*f0*t+phi)+5*sin(2*pi*f1*t+phi)+randn(1,N);
    [X(k,:),w]=periodogram(x(k,:),[],2^12);
end
X_ave=mean(X);

% Overlay of the periodograms
X_dB=20*log10(X');
figure;plot(w/(2*pi),X_dB)
axis([0 1/2 -80 10+max(max(X_dB))])

% Average of the periodograms
X_ave_dB=20*log10(X_ave);
figure;plot(w/(2*pi),X_ave_dB)
axis([0 1/2 -80 10+max(X_ave_dB)])

%% Periodogram of the white noise
clc
clear all;

fs=25e2; % Sampling frequency
Ts=1/fs; % Samplig period
f0=0.2*fs; % The frequency of the first sinusoid
f1=0.21*fs; % The frequency of the second sinusoid

N=256; % Number of samples in the signal
t=0:Ts:(N-1)*Ts; % The time vector


% Generating 50 sample signals and their periodograms
for k=1:50
    x(k,:)=randn(1,N);
    [X(k,:),w]=periodogram(x(k,:),[],2^12);
end
X_ave=mean(X);

% Overlay of the periodograms
X_dB=20*log10(X');
figure;plot(w/(2*pi),X_dB)
axis([0 1/2 -80 10+max(max(X_dB))])

% Average of the periodograms
X_ave_dB=20*log10(X_ave);
figure;plot(w/(2*pi),X_ave_dB)
axis([0 1/2 -80 10+max(X_ave_dB)])
