clc, clf, clear, close all
load('done_part3.mat')

%% General
N = 8192;
dt = t(2)-t(1);
T = t(end)+dt;

f0 = 1/T;
fs = 1/dt;

s = tf('s');
w = freq*2*pi;

opts = optimset('Display','off');

%% Part 4 - Extended Model
omega = ft_pkslocs*f0;
Hpe = magHp;

%% Model D

% do the exact same as in part3.m but now for the extended model

% x = [Kp, Ti, tau, Wnm, Znm]
x_d = zeros(5,5);
x_d(1,:) = [ 2    0.8 0.35   1.4 0.5];
x_d(2,:) = [ 3    1.6 1.25   3   0.2];
x_d(3,:) = [ 2.5  4.7 1.58   2   0.4];
x_d(4,:) = [ 2.75 3.2 2.83   0.5 0.5];
x_d(5,:) = [ 2.25 2.4 0.72   4   0.3];

lb       = [0.01 0.1 0.01  0.1  0.1];
ub       = [2.5  4   4     2    0.75];

Pars_D = zeros(5,5);
Cost_D = zeros(5,1);


J_D = @(x) sum(abs(Hpe-x(1)./(x(2).*1i.*omega+1).*exp(-1i.*omega.*x(3)).*(x(4).^2)./((1i.*omega).^2+2.*x(4).*x(5).*1i.*omega+x(4).^2)).^2);

for j=1:5
    [Pars_D(j,:),Cost_D(j,:)] = fmincon(J_D,x_d(j,:),[],[],[],[],lb,ub,[],opts);
end

idxD = find(Cost_D==min(Cost_D));
D_opt = Pars_D(idxD,:);

sys_D = D_opt(1)/(D_opt(2)*s+1)*exp(-s*D_opt(3))*(D_opt(4)^2)/(s^2+2*D_opt(4)*D_opt(5)*s+D_opt(4)^2);
[magD, phaD, wD] = bode(sys_D);


% Find system response
Kp_D    = D_opt(1);
Ti_D    = D_opt(2);
tau_D   = D_opt(3);
Wnm_D   = D_opt(4);
Znm_D   = D_opt(5);

res_D = zeros(length(w),1);
for k=1:length(w)
    res_D(k,1) = Kp_D  *  exp(-1i*w(k)*tau_D)  /  (1i*w(k)*Ti_D+1)  *  (Wnm_D^2) / ( (1i*w(k))^2 + 2*(1i*w(k))*Wnm_D*Znm_D + Wnm_D^2 );
end
mag_D = abs(res_D);
pha_D = rad2deg(unwrap(angle(res_D)));

%% Plots
figure(1)
sgtitle('Parameter Estimation');
subplot(2,1,1);
loglog(omega,magHp,'ko'); hold on
loglog(w,mag_A,'b');
loglog(w,mag_B,'g');
loglog(w,mag_C,'r');
loglog(w,mag_D,'k');
legend('H_{p} Estimation', 'Model A', 'Model B', 'Model C', 'Model D','location','northwest');
xlim([omega(1)/2 omega(end)*2]); grid on; ylim([8e-2 1e3]);
xlabel('\omega (rad/s)'); ylabel('|H(j\omega)| (abs)')

subplot(2,1,2)
semilogx(omega,phaHp,'ko'); hold on
semilogx(w,pha_A,'b');
semilogx(w,pha_B,'g');
semilogx(w,pha_C,'r');
semilogx(w,pha_D,'k');
legend('H_{p} Estimation', 'Model A', 'Model B', 'Model C', 'Model D','location','southwest');
xlim([omega(1)/2 omega(end)*2]); grid on; ylim([-600 100]);
xlabel('\omega (rad/s)'); ylabel('\angle H(j\omega) (deg)')

%% Print results
fprintf('Optimal sets found: \n\n')
fprintf('    Kp    tau    Wnm    Znm     Tl    Ti \n');
fprintf('A %-.4f %-.4f %-.4f %-.4f \n', Kp_A, tau_A, Wnm_A, Znm_A);
fprintf('B %-.4f %-.4f %-.4f %-.4f %-.4f \n', Kp_B, tau_B, Wnm_B, Znm_B,Tl_B);
fprintf('C %-.4f %-.4f %-.4f %-.4f %-.4f %-.4f \n', Kp_C, tau_C, Wnm_C, Znm_C, Tl_C, Ti_C);
fprintf('D %-.4f %-.4f %-.4f %-.4f        %-.4f \n', Kp_D, tau_D, Wnm_D, Znm_D, Ti_D);

save('done_part4.mat')