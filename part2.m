clc, clf, clear, close all
load('ae4320_dataset.mat')

%%% General
dt = t(2)-t(1);         % time step
T = t(end)+dt;          % sampling time
N = length(t);          % number of samples

f0 = 1/T;               % frequency resolution
fs = 1/dt;              % sampling frequency
fmax = N/2*f0;          % max frequency

%%% Part 2 - Identifying the pilot control dynamics FRF
%% Determine Excitation Frequencies of Forcing Function
% The time signals are plotted
% The double sideband fourier transform is plotted
% The single sideband fourier transform is plotted
%
% The excitation frequencies can be determined by looking at the peaks of
% the fft

X = fft(ft); U = abs(X/N); L = U(1:N/2+1); L(2:end-1) = 2*L(2:end-1);
Y = fft(fd); V = abs(Y/N); M = V(1:N/2+1); M(2:end-1) = 2*M(2:end-1);

% figure(1)
% sgtitle('Time Signals f_{t}(t) and f_{d}(t)')
% subplot(1,2,1); plot(t,ft); title('f_{t}(t)'); xlabel('t (s)');
% subplot(1,2,2); plot(t,fd); title('f_{d}(t)'); xlabel('t (s)');
% 
% figure(2)
f = 0:f0:2*fmax-f0;
% sgtitle('FFT (dsb)')
% subplot(1,2,1); plot(f,U); title('fft(f_{t})'); xlabel('f (Hz)');
% subplot(1,2,2); plot(f,V); title('fft(f_{d})'); xlabel('f (Hz)');
% 
% figure(3)
f = 0:f0:fmax;
% sgtitle('FFT (ssb)')
% subplot(1,2,1); plot(f,L); title('fft(f_{t})'); xlabel('f (Hz)'); xlim([0 50]);
% subplot(1,2,2); plot(f,M); title('fft(f_{d})'); xlabel('f (Hz)'); xlim([0 50]);

% Peak frequencies 
% (This very neat function requires Signal Processing Toolbox)
[pks,loc] = findpeaks(L,f/f0,'NPeaks',10,'Threshold',1e-2); 
f_peak_ft = loc*f0; % peak frequencies of f_t (Hz)
a_peak_ft = pks;    % peak amplitudes  of f_t (--)

[pks,loc] = findpeaks(M,f/f0,'NPeaks',10,'Threshold',1e-2);
f_peak_fd = loc*f0; % peak frequencies of f_d (Hz)
a_peak_fd = pks;    % peak amplitudes  of f_d (--)

f_peak_all = sort([f_peak_ft, f_peak_fd]/f0); % combined peak frequencies (idx)
f_peak_min = min(f_peak_all);
f_peak_max = max(f_peak_all);

%% Estimate Frequency Response Function of the Pilot control Dynamics
%
%  ^         U(jw)
%  H_p(jw) = -----
%            E(jw) 
%
% This is only done at the excitation frequencies of U and E
% So first a fourier analysis needs to be done again :)

e = identification.e; u = identification.u;
E = fft(e); F = abs(E/N); G = F(1:N/2+1); G(2:end-1) = 2*G(2:end-1);
U = fft(U); V = abs(U/N); W = V(1:N/2+1); W(2:end-1) = 2*W(2:end-1);

% figure(4)
% sgtitle('Time Signals e(t) and u(t)')
% subplot(1,2,1); plot(t,e); title('e(t)'); xlabel('t (s)');
% subplot(1,2,2); plot(t,u); title('u(t)'); xlabel('t (s)');

% figure(5)
f = 0:f0:2*fmax-f0;
% sgtitle('FFT (dsb)')
% subplot(1,2,1); plot(f,F); title('fft(e)'); xlabel('f (Hz)');
% subplot(1,2,2); plot(f,V); title('fft(u)'); xlabel('f (Hz)');
% 
% figure(6)
f = 0:f0:fmax;
% sgtitle('FFT (ssb)')
% subplot(1,2,1); plot(f,G); title('fft(e)'); xlabel('f (Hz)');
% subplot(1,2,2); plot(f,W); title('fft(u)'); xlabel('f (Hz)');


% magE = G(f_peak_all);
% phaE = angle(E/N)*180/pi;
% phaE = phaE(1:N/2+1);
% phaE = phaE(f_peak_all);
%  
% magU = W(f_peak_all);
% phaU = angle(U/N)*180/pi;
% phaU = phaU(1:N/2+1);
% phaU = phaU(f_peak_all);
% 
% magHp = magU/magE;
% phaHp = phaU-phaE;

Seft = (E(1:N)).*(X(1:N))/N;
Suft = (U(1:N)).*(Y(1:N))/N;

omegaf = f_peak_all*f0*2*pi;
Hp = Suft(f_peak_all)./Seft(f_peak_all);
magHp = abs(Hp);
phaHp = rad2deg(unwrap(angle(Hp)));

fprintf('This is an estmiation of the pilots control behaviour. \n')
fprintf('It only gives a real estimation at the frequencies of the peaks in the second ssb graph. \n')
figure(7)
sgtitle('H_p(jw)')
subplot(2,1,1); loglog(2*pi*f_peak_all*f0, magHp,'Color','#0072BD'); ylim([5e-5 5e-3]); yticks([1e-5 1e-4 1e-3]);
title('Magnitude'); xlabel('Frequency (rad/s)'); ylabel('Magnitude (abs)');
subplot(2,1,2); semilogx(2*pi*f_peak_all*f0, phaHp); ylim([-200 +200]);
title('Phase'); xlabel('Frequency (rad/s)'); ylabel('Phase (deg)');

save('part2_done.mat')