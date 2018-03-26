function [Px,w]=music(x,p,M,N)
% MUSIC Algorithm for Frequency Estimation
% By Ashkan Ashrafi
% San Diego, CA
% May 2nd, 2014
%
% INPUTS:
% x         :      Input signal
% p         :      Number of frequencies
% M         :      Number of available samples of the autocorrelation
% N         :      Size of frequency grid. It has be a power of 2
% 
% OUTPUTS:
% P         :      Frequency estimation function

x=x(:);
if M<p+1 || length(x)<M
    error('Size of R is not correct');
end
Lx=length(x);
x=x-mean(x);
rx=1/Lx*xcorr(x,x);
R=toeplitz(rx(Lx:Lx+M-1));
[V,d]=eig(R);
[ds,Ix]=sort(diag(d),'ascend');
Pxx=abs(fft(V(:,Ix(1:M-p)),N));
Px=-10*log10(sum(Pxx,2));
Px=Px(2:N/2);
w=2*pi/N:2*pi/N:(N/2-1)*2*pi/N;























