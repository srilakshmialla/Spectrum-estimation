function [Px,w]=phdm(x,p,N)
% Function for Pisarenko Harmonic Decompostion
%
% INPUTS:
% x         :      Input signal
% p         :      Number of frequencies
% N         :      Size of frequency grid
% 
% OUTPUTS:
% P         :      Frequency estimation function


x=x(:);
M=length(x);
x=x-mean(x);
rx=1/M*xcorr(x,x);
R=toeplitz(rx(M:M+p));
[V,d]=eig(R);
[ds,Ix]=sort(diag(d),'ascend');
Px=-20*log10(abs(fft(V(:,Ix(1)),N)));
Px=Px(2:N/2);
w=2*pi/N:2*pi/N:(N/2-1)*2*pi/N;