function [power,confd,freq,time] = pwavelet(x,fs,harm,wave,sigc)
% TWAVELET - proposed wavelet transform

% 1. pre-process data
sigma2 = std(x)^2; dt = 1/fs;
x = (x - mean(x))/sqrt(sigma2);
time = (0:dt:length(x)*dt-dt)';

% 2. construct wavenumber array
n = length(x);
wk = (1:n/2).*((2.*pi*fs/n));
wk = [0., wk, -wk(fix((n-1)/2):-1:1)];      % [Eqn(5)]
X = fft(x,length(x));                       % [Eqn(3)]

% 3. scale design
ds = harm.df/10;
s0 = 1/harm.fh;                                  % smallest possible
s1 = log2(1/(s0*harm.fl));
s = s0 * 2.^(0:ds:s1);                      % [Eqn(9)]

% 4. wavelet transform
Hw = (wk > 0.);                             % Heaviside step function
switch(wave.basis)                               %[Table(1)]
    case 'morlet'
        wk0 = 6; dof = 2;
        lambda = (4*pi)/(wk0 + sqrt(2 + wk0^2));
        for i=1:length(s)
            hpsi_0 = (pi^(-0.25)).*Hw.*exp(-(s(i).*wk - wk0).^2/2);
            hpsi_n(:,i) = sqrt(s(i)*2*pi*fs) .* hpsi_0;           %[Eqn(6)]
            Wn(i,:) = ifft(X.*hpsi_n(:,i));                       %[Eqn(4)]
        end
        
    case 'paul'
        m = 4; dof = 2;
        lambda = 4*pi/(2*m+1);
        for i=1:length(s)
            hpsi_0 = 2^m/sqrt(m*prod(2:(2*m-1))).*Hw.*((s(i).*wk).^m).*exp(-s(i).*wk);
            hpsi_n(:,i) = sqrt(s(i)*2*pi*fs) .* hpsi_0;           %[Eqn(6)]
            Wn(i,:) = ifft(X.*hpsi_n(:,i));                       %[Eqn(4)]
        end
        
    case 'dog'
        m = 2; dof = 1;
        lambda = 2*pi/sqrt(m+1/2);
        for i=1:length(s)
            hpsi_0 = 1i^m/sqrt(gamma(m+0.5)).*((s(i).*wk).^m).*exp(-(s(i).*wk).^2/2);
            hpsi_n(:,i) = sqrt(s(i)*2*pi*fs) .* hpsi_0;           %[Eqn(6)]
            Wn(i,:) = ifft(X.*hpsi_n(:,i));                       %[Eqn(4)]
        end
end
power.c = (abs(Wn)).^2 ;                                  % power spectrum
freq = 1./(lambda*s)';

% 4. Significance Level
Pk = (1-sigc.alpha^2) ./ ...
     (1+sigc.alpha^2-2*sigc.alpha*cos(freq*2*pi));         % [Eqn(16)]
chi2 = chi2inv(sigc.level,dof);
sig = Pk*chi2/dof;                                    % [Eqn(18)]
confd.c = power.c ./ (sig*ones(1,n));                  % Confidence level
% ratio > 1 = significant power

% 5. Cone-of-influence
coi = lambda/sqrt(2);             % cone-of-influence [Sec.3g]
coi = fs./(coi*[1E-1,1:((n+1)/2-1),fliplr((1:(n/2-1))),1E-1]);

%% Significance
switch(wave.basis)                               %[Table(2)]
    case 'morlet'
        Cdelta = 0.776;
        Gamma = 2.32;
        ds0 = 0.60;
    case 'paul'
        Cdelta = 1.132;
        Gamma = 1.17;
        ds0 = 1.5;
    case 'dog'
        Cdelta = 3.541;
        Gamma = 1.43;
        ds0 = 1.4;
end

%% 6. time-averaging
% 1. Average Power
power.f = sigma2*(sum(power.c')/n);

% 2. Significance level
Pk = sigma2*Pk;
dofn = n - s;
truncate = find(dofn < 1);
dofn(truncate) = ones(size(truncate));
dofn = dof*sqrt(1 + (dofn*dt/Gamma ./ s).^2 );   % [Eqn(23)]
truncate = find(dofn < dof);
dofn(truncate) = dof*ones(size(truncate));
for i = 1:length(dofn)
    chi2(i) = chi2inv(sigc.level,dofn(i));
    sig(i) = Pk(i)*chi2(i)/dofn(i);
end
confd.f = (mean(sig));% - min(power.f))./std(power.f);
%power.f = (power.f - min(power.f))./std(power.f);

%% Scale-Averaged
% 1. Averaged Power
s_n = (s'*ones(1,n));
a_idx = find((s >= 1/sigc.range(2)) & (s < 1/sigc.range(1)));
power.t = sigma2*ds*dt/Cdelta * sum(power.c(a_idx,:)./s_n(a_idx,:));          % [Eqn(24)]

% 2. Significance Level
na = length(a_idx);
Savg = 1./sum(1 ./ s(a_idx));                           % [Eqn(25)]
Smid = exp((log(1/sigc.range(2))+log(sigc.range(1)))/2);                % power-of-two midpoint
dofn = (dof*na*Savg/Smid)*sqrt(1 + (na*ds/ds0)^2);    % [Eqn(28)]
Pk = Savg*sum(Pk(a_idx) ./ s(a_idx)');                   % [Eqn(27)]
chi2 = chi2inv(sigc.level,dofn);
confd.t = (ds*dt/Cdelta/Savg)*Pk*chi2/dofn;               % [Eqn(26)]

%% Wavelet Reconstruction
for i=1:length(s)
    Wd(i) = (1/length(Wn))*sum(hpsi_n(:,i));    %[Eqn(12)]
end
Cd = sum(real(Wd)./sqrt(s));                    %[Eqn(13)]
sn = repmat(s(:),[1,size(Wn,2)]);
xr = (1/Cd)*sum(real(Wn)./sqrt(sn),1);          %[Eqn(11)]


%% Wavelet results overview
if wave.figs == 1
    figure('units','normalized','outerposition',[0 0 1 1])
    hax1 = subplot('position',[0.08 0.77 0.65 0.2]);
    plot(time,x)
        xlim([time(1),time(end)])
        set(gca,'XtickLabel',[])
        ylabel('Data')

    hax2 = subplot('position',[0.08 0.34 0.65 0.4]);
    contour(time,freq,power.c,25); set(gca,'yscale','log')
        xlim([time(1),time(end)])
        ylim([min(freq),max(freq)])
        set(gca,'XtickLabel',[])
        ylabel('Freq'), hold on
    contour(time,freq,confd.c,[-1e4,1],'k','Linewidth',1.25);
    semilogy(time,coi,'k')

    hax3 = subplot('position',[0.75 0.34 0.2 0.4]);
    semilogy(power.f,freq), hold on
    semilogy(ones(length(freq),1)*confd.f,freq,'--k','Linewidth',1.25)
        xlim([0,1.25*max(power.f)])
        ylim([min(freq),max(freq)])
        set(gca,'YtickLabel',[]), xlabel('Normalized Power')
    
    hax4 = subplot('position',[0.08 0.11 0.65 0.2]);
    plot(time,power.t)
        xlim([time(1),time(end)]), xlabel('Time')
        ylabel('Power'), hold on
    plot([time(1),time(end)],[confd.t,confd.t],'--k','Linewidth',1.25)

    linkaxes([hax2,hax3], 'y');
    linkaxes([hax1,hax2,hax4], 'x');
end

end

