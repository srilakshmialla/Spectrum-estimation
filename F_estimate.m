% Advanced EE=658
% Frequency estimation 
% By : Ashkan Ashrafi
% San Diego, CA
% 5/1/2014

close all;
clear all;
clc;

fs=25e3; % Sampling frequency
Ts=1/fs; % Samplig period
f0=800;
f1=2300;
f2=5900;
f3=8700;



Nf=1024;
p=4;

%% Signal
N=512;
t=0:Ts:(N-1)*Ts; % The time vector
A0=2.1;
A1=1.9;
A2=2.6;
A3=1.1;
An=5;
for k=1:10
    phi1= -pi + 2*pi.*rand(1,1);
    phi2= -pi + 2*pi.*rand(1,1);
    phi3= -pi + 2*pi.*rand(1,1);
    phi4= -pi + 2*pi.*rand(1,1);
    x(k,:)=A0*sin(2*pi*f0*t+phi1)+A1*sin(2*pi*f1*t+phi2)+A2*sin(2*pi*f2*t+phi3)...
        +A3*sin(2*pi*f3*t+phi4)+An*randn(1,N);
end
SNR=10*log10(A0^2+A1^2+A2^2+A3^2)/An^2

%% Minimum Variance

for k=1:10
    [X(k,:),wmv]=MV(x(k,:),p,Nf);
    Xmv(k,:)=10*log10(abs(X(k,:)));
end

%% Blackman-Tukey

M=Nf/4;
wind=hamming(2*M-1);
for k=1:10
    [X(k,:),wbt]=per_smooth(x(k,:),wind,M,Nf);
    Xbt=10*log10(abs(X(k,:)));
end


%% Pisarenko

for k=1:10
    [Px_p(k,:),wp]=phdm(x(k,:),p,Nf);
end

%% MUSIC

M=64;
for k=1:10
    [Px_m(k,:),wm]=music(x(k,:),p,M,Nf);
end

%% Miniumum Norm

M=64;
for k=1:10
    [Px_mn(k,:),wmn]=min_norm(x(k,:),p,M,Nf);
end

%% Plots

% close all;

figure;plot(t,x(1,:))
xlabel('Time (Sec)')
title('The Signal')


figure;plot(fs*wbt/(2*pi),Xmv')
xlabel('Frequency (Hz)')
ylabel('Px');
title('Minimum Variance Method')

figure;plot(fs*wmv/(2*pi),Xbt')
xlabel('Frequency (Hz)')
ylabel('Px');
title('Blackma-Tukey Method')

figure;plot(fs*wp/(2*pi),Px_p')
xlabel('Frequency (Hz)')
ylabel('Px');
title('Pisarenko Method')

figure;plot(fs*wm/(2*pi),Px_m')
xlabel('Frequency (Hz)')
ylabel('Px');
title('MUSIC Method')

figure;plot(fs*wmn/(2*pi),Px_mn')
xlabel('Frequency (Hz)')
ylabel('Px');
title('Minimum Norm Method')
