clc, clf, clear, close all
load('part3_done.mat')

%% General
N = 8192;
dt = t(2)-t(1);
T = t(end)+dt;

f0 = 1/T;
fs = 1/dt;

s = tf('s');

% magHp phaHp

%% Part 4 - Extended Model
omega = ft_pkslocs*f0;
Hpe = magHp;

%% Model D

% do the exact same as in part3.m but now for the extended model

% x = [Kp, Ti, tau, Wnm, Znm]
x_d = zeros(5,5);
x_d(1,:) = [0, 1, 1, 1, 1];
x_d(2,:) = [1, 0, 1, 1, 1];
x_d(3,:) = [1, 1, 0, 1, 1];
x_d(4,:) = [1, 1, 1, 0, 1];
x_d(5,:) = [1, 1, 1, 1, 0];

Pars_D = zeros(5,5);
Cost_D = zeros(5,1);

for k=1:10
    J_D = @(x) sum(abs(Hpe(k)-x(1)/(x(2)*1i*omega(k)+1)*exp(-1i*omega(k)*x(3))*(x(4)^2)/((1i*omega(k)^2)+2*x(4)*x(5)*1i*omega(k)+x(4)^2))^2);
end

for j=1:5
    Pars_D(j,:) = fminsearch(J_D,x_d(j,:));
    Cost_D(j,1) = J_D(Pars_D(j,:));
end

idxD = find(Cost_D==min(Cost_D));
D_opt = Pars_D(idxD,:);

sys_D = D_opt(1)/(D_opt(2)*s+1)*exp(-s*D_opt(3))*(D_opt(4)^2)/(s^2+2*D_opt(4)*D_opt(5)*s+D_opt(4)^2);
[magD, phaD, wD] = bode(sys_D);

%% Plots
figure(1)
subplot(2,1,1);
loglog(omega,magHp,'k'); hold on
loglog(wA,reshape(magA,[length(magA),1]),'blue--'); 
loglog(wB,reshape(magB,[length(magB),1]),'green--'); 
loglog(wC,reshape(magC,[length(magC),1]),'red--');
loglog(wD,reshape(magD,[length(magD),1]),'yellow--');
legend('H_{p} Estimation', 'Model A', 'Model B', 'Model C', 'Model D');
xlim([omega(1) omega(end)]);

subplot(2,1,2)
semilogx(omega,phaHp,'k'); hold on
semilogx(wA,reshape(phaA,[length(phaA),1]),'blue--'); 
semilogx(wB,reshape(phaB,[length(phaB),1]),'green--'); 
semilogx(wC,reshape(phaC,[length(phaC),1]),'red--');
semilogx(wD,reshape(phaD,[length(phaD),1]),'yellow--');
legend('H_{p} Estimation', 'Model A', 'Model B', 'Model C', 'Model D');
xlim([omega(1) omega(end)]);

save('part4_done.mat')