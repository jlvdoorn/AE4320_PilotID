clc, clf, clear, close all
load('part4_done.mat')

%% General
N = 8192;
dt = t(2)-t(1);
T = t(end)+dt;

f0 = 1/T;
fs = 1/dt;

freq  = f0*(1:N/2);
omega = freq*2*pi;

ev = validation.e;
uv = validation.u;

%% Part 5 - Model Validation
ev_dft  = fft(ev,N);
uv_dft  = fft(uv,N);

Hpv = u_dft(ft_pkslocs)./e_dft(ft_pkslocs);
magHpv = abs(Hpv); phaHpv = rad2deg(unwrap(angle(Hpv)));

figure(1)
sgtitle('Validation of Hp(j\omega)')

subplot(2,1,1); loglog(ft_pkslocs*f0, magHpv,'ob'); hold on;
loglog(ft_pkslocs*f0, magHp,'*g');
xlim([5e-2 5e0]); ylim([1e0 1e2]); grid on;
legend('Real H_{p}(j\omega)','Estimated H_{p}(j\omega)');
xlabel('\omega (rad/s)'); ylabel('|Hp(j\omega)| (abs)')

subplot(2,1,2); semilogx(ft_pkslocs*f0,phaHpv,'ob'); hold on;
semilogx(ft_pkslocs*f0, phaHp,'*g');
xlim([5e-2 5e0]); ylim([-450 90]); grid on;
xlabel('\omega (rad/s)'); ylabel('\angle Hp(j\omega) (deg)')
legend('Real H_{p}(j\omega)','Estimated H_{p}(j\omega)');