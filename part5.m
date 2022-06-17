clc, clf, clear, close all
load('done_part4.mat')

%% General
N = 8192;
dt = t(2)-t(1);
T = t(end)+dt;

f0 = 1/T;
fs = 1/dt;

freq  = f0*(1:N/2);
omega = freq*2*pi;

%% Part 5 - Model Validation
ev = validation.e;
uv = validation.u;

%% FRF of validation set
ev_dft = fft(ev,N);
uv_dft = fft(uv,N);
ft_dft = fft(ft,N);

[ft_pks, ft_pkslocs] = findpeaks(abs(ft_dft(1:N/2)),'MinPeakHeight',1e-1);

% Hp=U(jw)/E(jw) but now based on validation dataset
Hpv = uv_dft(ft_pkslocs)./ev_dft(ft_pkslocs);
magHpv = abs(Hpv); phaHpv = rad2deg(unwrap(angle(Hpv)));

figure(1)
sgtitle('Validation of Hp(j\omega)')

subplot(2,1,1); loglog(ft_pkslocs*f0, magHpv,'*g'); hold on;
loglog(ft_pkslocs*f0, magHp,'ob');
loglog(w,mag_A,'b');
loglog(w,mag_B,'g');
loglog(w,mag_C,'r');
loglog(w,mag_D,'k');
xlim([5e-2 5e0]); ylim([1e0 1e2]); grid on; ylim([1e-1 1e2]);
xlabel('\omega (rad/s)'); ylabel('|Hp(j\omega)| (abs)')
legend('Real H_{p}(j\omega)','Estimated H_{p}(j\omega)','Model A','Model B','Model C','Model D','NumColumns',3,'location','southwest');

subplot(2,1,2); semilogx(ft_pkslocs*f0,phaHpv,'*g'); hold on;
semilogx(ft_pkslocs*f0, phaHp,'ob');
semilogx(w,pha_A,'b');
semilogx(w,pha_B,'g');
semilogx(w,pha_C,'r');
semilogx(w,pha_D,'k');
xlim([5e-2 5e0]); ylim([-450 90]); grid on; ylim([-600 100]);
xlabel('\omega (rad/s)'); ylabel('\angle Hp(j\omega) (deg)')
legend('Real H_{p}(j\omega)','Estimated H_{p}(j\omega)','Model A','Model B','Model C','Model D','NumColumns',3,'location','southwest');

%% Cost Function of Model A-D based on real Hp
omega = ft_pkslocs*f0;
Hpv = magHpv;

% create cost functions over all peak frequencies
J_Av = @(x) sum(abs(Hpv-x(1).* exp(-1i*omega*x(2)).*(x(3).^2)./((1i*omega).^2 + 2.*x(3).*x(4).*1i.*omega + x(3).^2)).^2);     
J_Bv = @(x) sum(abs(Hpv-x(1).*(x(2).*1i.*omega+1).*exp(-1i.*omega.*x(3)).*(x(4).^2)./((1i.*omega).^2+2.*x(4).*x(5).*1i.*omega+x(4).^2)).^2);
J_Cv = @(x) sum(abs(Hpv-x(1).*(x(2).*1i.*omega+1)./(x(3).*1i.*omega+1).*exp(-1i.*omega.*x(4)).*(x(5).^2)./((1i.*omega).^2+2.*x(5).*x(6).*1i.*omega+x(5).^2)).^2);
J_Dv = @(x) sum(abs(Hpv-x(1)./(x(2).*1i.*omega+1).*exp(-1i.*omega.*x(3)).*(x(4).^2)./((1i.*omega).^2+2.*x(4).*x(5).*1i.*omega+x(4).^2)).^2);


% calculate cost for all cost functions based on the opt. parameter set
Cost_Av = J_Av(A_opt);
Cost_Bv = J_Bv(B_opt);
Cost_Cv = J_Cv(C_opt);
Cost_Dv = J_Dv(D_opt);

save('done_part5.mat')