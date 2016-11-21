function [nt,nf,power,ta,fa] = pwavelet_toyota(time,xn,fr,df,dt,fa)
% Author: Saito, Toyota company
% Notes:
% it seems to me that this is based on local fft insteaed of wavelet
% wavelet should not have a fixed DT (see typical fft vs. wavelet figures)
% the convolution function used here seems to be sin/cos (not finite)
% the data is first fourier transformed (xn.*exp((t-t1)*2*pi*f))
% and then multiplied in frequency domain with cos/sin functions

fl = fr(1); fh = fr(2);
Ts = time(2);
fs = 1/Ts;
Te = time(end);
nrofs = length(xn);
nroff = floor((fh - fl)/df);
nroft = floor(Te/dt+1);

tauSet = 5;
switch tauSet
    case 2, tauN = 8; unitV = 28.2*10;
    case 3, tauN = 16; unitV = 39.45*10;
    case 4, tauN = 32; unitV = 56.34*10;
    case 5, tauN = 64; unitV = 79.75*10;
end

mwave1 = zeros(nrofs,nroft);
mwave2 = zeros(nrofs,nroft);
power = zeros(nroft,nroff);
ma = zeros(nroff,1);

for m = 1:(fh - fl)/df;
    f = fl + (m-1)*df;
    for n = 1:Te/dt+1;
        tt = (time - (n-1)*dt)*2*pi*f;
        mwave1(:,n) = xn.* ( exp(-(tt).^2/tauN).*cos(tt) );
        mwave2(:,n) = xn.* ( exp(-(tt).^2/tauN).*sin(tt) );
        power(n,m) = sqrt((sum(mwave1(:,n)).^2+sum(mwave2(:,n)).^2)/2*f^2)/unitV*fs;
        nt(:,n) = (n - 1);
    end
    nf(:,m) = f;
end
nt = nt*dt;

for i=1:length(nf)
    ta(i) = mean(power(:,i));
end

idx0 = find(nf > fa(1),1);
idx1 = find(nf > fa(2),1);
for i=1:length(nt)
    fa(i) = mean(power(i,idx0:idx1));
end

end
