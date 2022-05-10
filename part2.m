clc, clf, clear, close all
load('ae4320_dataset.mat')

%% General
N = 8192;
dt = t(2)-t(1);
T = t(end)+dt;

f0 = 1/T;
fs = 1/dt;

freq  = f0*(1:N/2);
omega = freq*2*pi;

e = identification.e;
u = identification.u;

%% Part 2 - Identifying the pilot control dynamics FRF
ft_dft = fft(ft,N);
fd_dft = fft(fd,N);

[ft_pks, ft_pkslocs] = findpeaks(abs(ft_dft(1:N/2)),'MinPeakHeight',1e-1);
[fd_pks, fd_pkslocs] = findpeaks(abs(fd_dft(1:N/2)),'MinPeakHeight',1e-1);

figure(1) % dft(ft) and dft(fd)
sgtitle('DFT of the Forcing Functions')
subplot(1,2,1);
loglog(omega, abs(ft_dft(1:N/2)),'b'); title('ft'); hold on
loglog(omega(ft_pkslocs), ft_pks,'ob'); xlim([7e-2 3e2]); ylim([1e-6 1e4]);
xlabel('\omega (rad/s)'); ylabel('|dft(ft) (abs)|');

subplot(1,2,2);
loglog(omega, abs(fd_dft(1:N/2)),'b'); title('fd'); hold on
loglog(omega(fd_pkslocs), fd_pks,'ob'); xlim([7e-2 3e2]); ylim([1e-6 1e4]);
xlabel('\omega (rad/s)'); ylabel('|dft(fd) (abs)|');

e_dft  = fft(e,N);
u_dft  = fft(u,N);

% figure(2) % dft(e) and dft(u)
% loglog(omega, abs(e_dft(1:N/2)),'b'); hold on
% loglog(omega(fd_pkslocs), abs(e_dft(fd_pkslocs)),'ob');
% loglog(omega(ft_pkslocs), abs(e_dft(ft_pkslocs)),'ob');
% 
% loglog(omega, abs(u_dft(1:N/2)),'g'); 
% loglog(omega(fd_pkslocs), abs(u_dft(fd_pkslocs)),'og');
% loglog(omega(ft_pkslocs), abs(u_dft(ft_pkslocs)),'og');
% 
% Seft = (e_dft(1:N)).*(ft_dft(1:N))/N;
% Suft = (u_dft(1:N)).*(ft_dft(1:N))/N;
% 
% figure(3) % Seft and Suft
% loglog(omega, abs(Seft(1:N/2)),'b'); hold on
% loglog(omega(ft_pkslocs), abs(Seft(ft_pkslocs)),'ob');
% 
% loglog(omega, abs(Suft(1:N/2)),'g'); hold on
% loglog(omega(ft_pkslocs), abs(Suft(ft_pkslocs)),'og');

Hp = u_dft(ft_pkslocs)./e_dft(ft_pkslocs);
magHp = abs(Hp); phaHp = rad2deg(unwrap(angle(Hp)));

figure(5) % Hp = U/E
sgtitle('Estimation of Hp(j\omega)')
subplot(2,1,1); loglog(ft_pkslocs*f0, magHp,'ob');
xlim([5e-2 5e0]); ylim([1e0 1e2]); grid on;
xlabel('\omega (rad/s)'); ylabel('|Hp(j\omega)| (abs)')
subplot(2,1,2); semilogx(ft_pkslocs*f0,phaHp,'ob');
xlim([5e-2 5e0]); ylim([-450 90]); grid on;
xlabel('\omega (rad/s)'); ylabel('\angle Hp(j\omega) (deg)')

save('part2_done.mat')