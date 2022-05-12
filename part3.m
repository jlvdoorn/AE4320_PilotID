clc, clf, clear, close all
load('part2_done.mat')

%% General
N = 8192;
dt = t(2)-t(1);
T = t(end)+dt;

f0 = 1/T;
fs = 1/dt;

s = tf('s');

% magHp phaHp

%% Part 3 - Parameter Estimation and Model Selection
omega = ft_pkslocs*f0;
Hpe = magHp;

%% Model A

% x = [Kp, tau, Wnm, Znm];
x_a = zeros(5,4);
x_a(1,:) = [0, 1, 1, 1];
x_a(2,:) = [1, 0, 1, 1];
x_a(3,:) = [1, 1, 0, 1];
x_a(4,:) = [1, 1, 1, 0];
x_a(5,:) = [1, 1, 1, 1];

Pars_A = zeros(5,4);
Cost_A = zeros(5,1);

for k=1:10
    J_A = @(x) sum(abs(Hpe(k)-x(1)*exp(-1i*omega(k)*x(2))*(x(3)^2)/((1i*omega(k)^2)+2*x(3)*x(4)*1i*omega(k)+x(3)^2))^2);
end

for j=1:5
    Pars_A(j,:) = fminsearch(J_A,x_a(j,:));
    Cost_A(j,1) = J_A(Pars_A(j,:));
end

idxA = find(Cost_A==min(Cost_A));
A_opt = Pars_A(idxA,:);

sys_A = A_opt(1)*exp(-s*A_opt(2))*(A_opt(3)^2)/(s^2+2*A_opt(3)*A_opt(4)*s+A_opt(3)^2);
[magA, phaA, wA] = bode(sys_A);

%% Model B

% x = [Kp, Tl, tau, Wnm, Znm];
x_b = zeros(5,5);
x_b(1,:) = [0, 1, 1, 1, 1];
x_b(2,:) = [1, 0, 1, 1, 1];
x_b(3,:) = [1, 1, 0, 1, 1];
x_b(4,:) = [1, 1, 1, 0, 1];
x_b(5,:) = [1, 1, 1, 1, 0];

Pars_B = zeros(5,5);
Cost_B = zeros(5,1);

for k=1:10
    J_B = @(x) sum(abs(Hpe(k)-x(1)*(x(2)*1i*omega(k)+1)*exp(-1i*omega(k)*x(3))*(x(4)^2)/((1i*omega(k)^2)+2*x(4)*x(5)*1i*omega(k)+x(4)^2))^2);
end

for j=1:5
    Pars_B(j,:) = fminsearch(J_B,x_b(j,:));
    Cost_B(j,1) = J_B(Pars_B(j,:));
end

idxB = find(Cost_B==min(Cost_B));
B_opt = Pars_B(idxB,:);

sys_B = B_opt(1)*(B_opt(2)*s+1)*exp(-s*B_opt(3))*(B_opt(4)^2)/(s^2+2*B_opt(4)*B_opt(5)*s+B_opt(4)^2);
[magB, phaB, wB] = bode(sys_B);

%% Model C

% x = [Kp, Tl, Ti, tau, Wnm, Znm];
x_c = zeros(5,6);
x_c(1,:) = [0, 1, 1, 1, 1, 1];
x_c(2,:) = [1, 0, 1, 1, 1, 1];
x_c(3,:) = [1, 1, 0, 1, 1, 1];
x_c(4,:) = [1, 1, 1, 0, 1, 1];
x_c(5,:) = [1, 1, 1, 1, 0, 1];

Pars_C = zeros(5,6);
Cost_C = zeros(5,1);

for k=1:10
    J_C = @(x) sum(abs(Hpe(k)-x(1)*(x(2)*1i*omega(k)+1)/(x(3)*1i*omega(k)+1)*exp(-1i*omega(k)*x(4))*(x(5)^2)/((1i*omega(k)^2)+2*x(5)*x(6)*1i*omega(k)+x(5)^2))^2);
end

for j=1:5
    Pars_C(j,:) = fminsearch(J_C,x_c(j,:));
    Cost_C(j,1) = J_C(Pars_C(j,:));
end

idxC = find(Cost_C==min(Cost_C));
C_opt = Pars_C(idxC,:);

sys_C = C_opt(1)*(C_opt(2)*s+1)/(C_opt(3)*s+1)*exp(-s*C_opt(4))*(C_opt(5)^2)/(s^2+2*C_opt(5)*C_opt(6)*s+C_opt(5)^2);
[magC, phaC, wC] = bode(sys_C);

%% Plots
figure(1)
subplot(2,1,1);
loglog(omega,magHp,'k'); hold on
loglog(wA,reshape(magA,[length(magA),1]),'blue--'); 
loglog(wB,reshape(magB,[length(magB),1]),'green--'); 
loglog(wC,reshape(magC,[length(magC),1]),'red--');
legend('H_{p} Estimation', 'Model A', 'Model B', 'Model C');
xlim([omega(1) omega(end)]);

subplot(2,1,2)
semilogx(omega,phaHp,'k'); hold on
semilogx(wA,reshape(phaA,[length(phaA),1]),'blue--'); 
semilogx(wB,reshape(phaB,[length(phaB),1]),'green--'); 
semilogx(wC,reshape(phaC,[length(phaC),1]),'red--');
legend('H_{p} Estimation', 'Model A', 'Model B', 'Model C');
xlim([omega(1) omega(end)]);

save('part3_done.mat')