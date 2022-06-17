clc, clf, clear, close all
load('done_part2.mat')

%% General
N = 8192;
dt = t(2)-t(1);
T = t(end)+dt;

f0 = 1/T;
fs = 1/dt;

s = tf('s');
w = freq*2*pi;

opts = optimset('Display','off');

%% Part 3 - Parameter Estimation and Model Selection
omega = ft_pkslocs*f0;
Hpe = Hp;

%% Model A

% define inco
% x = [Kp, tau, Wnm, Znm];
x_a = zeros(5,4);
x_a(1,:) = [ 2.5 0.5 2.1 0.11];
x_a(2,:) = [ 8.3 3.6 3.2 0.63];
x_a(3,:) = [ 9.5 5.8 9.6 0.75];
x_a(4,:) = [ 1.0 4.6 3.1 0.26];
x_a(5,:) = [ 2.5 3.5 1.5 0.10];

lb       = [0.1 0.01  0.1  0.1];
ub       = [10    10     10   0.75];

Pars_A = zeros(5,4);
Cost_A = zeros(5,1);

% create cost function over all peak frequencies 
J_A = @(x) sum(abs(  Hpe - x(1) .* exp(-1i*omega*x(2)) .* (x(3).^2) ./ ((1i*omega).^2 + 2.*x(3).*x(4).*1i.*omega + x(3).^2) ) .^2);     

% optimize cost function for all 5 incos
for j=1:5
    [Pars_A(j,:),Cost_A(j,:)] = fmincon(J_A,x_a(j,:),[],[],[],[],lb,ub,[],opts);
end

% find opt. parameter set based on minimum cost from 5 incos
idxA = find(Cost_A==min(Cost_A));
A_opt = Pars_A(idxA,:);

% define tf of model using opt. parameter set
sys_A = A_opt(1)*exp(-s*A_opt(2))*(A_opt(3)^2)/(s^2+2*A_opt(3)*A_opt(4)*s+A_opt(3)^2);
[magA, phaA, wA] = bode(sys_A, {5E-2 5E0});

% Find system response
Kp_A    = A_opt(1);
tau_A   = A_opt(2);
Wnm_A   = A_opt(3);
Znm_A   = A_opt(4);

res_A = zeros(length(w),1);
for k=1:length(w)
    res_A(k,1) = Kp_A  *  exp(-1i*w(k)*tau_A)  *  (Wnm_A^2) / ( (1i*w(k))^2 + 2*(1i*w(k))*Wnm_A*Znm_A + Wnm_A^2 );
end
mag_A = abs(res_A);
pha_A = rad2deg(unwrap(angle(res_A)));

%% Model B

% x = [Kp, Tl, tau, Wnm, Znm];
x_b = zeros(5,5);
x_b(1,:) = [ 2    0.8 0.35   1.4 0.5];
x_b(2,:) = [ 3    1.6 1.25   3   0.2];
x_b(3,:) = [ 2.5  4.7 1.58   2   0.4];
x_b(4,:) = [ 2.75 3.2 2.83   0.5 0.5];
x_b(5,:) = [ 2.25 2.4 0.72   4   0.3];

lb       = [0.01 0.1 0.01  0.1  0.1];
ub       = [2.5  5   3     2    0.75];

Pars_B = zeros(5,5);
Cost_B = zeros(5,1);

J_B = @(x) sum(abs(Hpe-x(1).*(x(2).*1i.*omega+1).*exp(-1i.*omega.*x(3)).*(x(4).^2)./((1i.*omega).^2+2.*x(4).*x(5).*1i.*omega+x(4).^2)).^2);

for j=1:5
    [Pars_B(j,:),Cost_B(j,:)] = fmincon(J_B,x_b(j,:),[],[],[],[],lb,ub,[],opts);
end

idxB = find(Cost_B==min(Cost_B));
B_opt = Pars_B(idxB,:);

sys_B = B_opt(1)*(B_opt(2)*s+1)*exp(-s*B_opt(3))*(B_opt(4)^2)/(s^2+2*B_opt(4)*B_opt(5)*s+B_opt(4)^2);
[magB, phaB, wB] = bode(sys_B, {5E-2 5E0});

% Find system response
Kp_B    = B_opt(1);
Tl_B    = B_opt(2);
tau_B   = B_opt(3);
Wnm_B   = B_opt(4);
Znm_B   = B_opt(5);

res_B = zeros(length(w),1);
for k=1:length(w)
    res_B(k,1) = Kp_B  *  exp(-1i*w(k)*tau_B)  *  (1i*w(k)*Tl_B+1)  *  (Wnm_B^2) / ( (1i*w(k))^2 + 2*(1i*w(k))*Wnm_B*Znm_B + Wnm_B^2 );
end
mag_B = abs(res_B);
pha_B = rad2deg(unwrap(angle(res_B)));

%% Model C

% x = [Kp, Tl, Ti, tau, Wnm, Znm];
x_c = zeros(5,6);
x_c(1,:) = [ 2    0.8 0.1 0.35   1.4 0.5];
x_c(2,:) = [ 3    1.6 0.5 1.25   3   0.2];
x_c(3,:) = [ 2.5  4.7 1.6 1.58   2   0.4];
x_c(4,:) = [ 2.75 3.2 0.8 2.83   0.5 0.1];
x_c(5,:) = [ 2.25 2.4 2.6 0.72   4   0.3];

lb       = [ 0.01 0.1 0.5 0.01  0.1  0.1];
ub       = [10    5   5   3    10    0.75];

Pars_C = zeros(5,6);
Cost_C = zeros(5,1);

J_C = @(x) sum(abs(Hpe-x(1).*(x(2).*1i.*omega+1)./(x(3).*1i.*omega+1).*exp(-1i.*omega.*x(4)).*(x(5).^2)./((1i.*omega).^2+2.*x(5).*x(6).*1i.*omega+x(5).^2)).^2);

for j=1:5
    [Pars_C(j,:),Cost_C(j,:)] = fmincon(J_C,x_c(j,:),[],[],[],[],lb,ub,[],opts);
end

idxC = find(Cost_C==min(Cost_C));
C_opt = Pars_C(idxC,:);

sys_C = C_opt(1)*(C_opt(2)*s+1)/(C_opt(3)*s+1)*exp(-s*C_opt(4))*(C_opt(5)^2)/(s^2+2*C_opt(5)*C_opt(6)*s+C_opt(5)^2);
[magC, phaC, wC] = bode(sys_C, {5E-2 5E0});

% Find system response
Kp_C    = C_opt(1);
Tl_C    = C_opt(2);
Ti_C    = C_opt(3);
tau_C   = C_opt(4);
Wnm_C   = C_opt(5);
Znm_C   = C_opt(6);

res_C = zeros(length(w),1);
for k=1:length(w)
    res_C(k,1) = Kp_C  *  exp(-1i*w(k)*tau_C)  *  (1i*w(k)*Tl_C+1)  /  (1i*w(k)*Ti_C+1)  *  (Wnm_C^2) / ( (1i*w(k))^2 + 2*(1i*w(k))*Wnm_C*Znm_C + Wnm_C^2 );
end
mag_C = abs(res_C);
pha_C = rad2deg(unwrap(angle(res_C)));

%% Plots
figure(1)
sgtitle('Parameter Estimation');
subplot(2,1,1);
loglog(omega,magHp,'ko'); hold on
loglog(w,mag_A,'b');
loglog(w,mag_B,'g');
loglog(w,mag_C,'r');
legend('H_{p} Estimation', 'Model A', 'Model B', 'Model C','location','northwest');
xlim([omega(1)/2 omega(end)*2]); grid on; ylim([8e-2 1e3]);
xlabel('\omega (rad/s)'); ylabel('|H(j\omega)| (abs)')

subplot(2,1,2)
semilogx(omega,phaHp,'ko'); hold on
semilogx(w,pha_A,'b');
semilogx(w,pha_B,'g');
semilogx(w,pha_C,'r');
legend('H_{p} Estimation', 'Model A', 'Model B', 'Model C','location','southwest');
xlim([omega(1)/2 omega(end)*2]); grid on; ylim([-600 100]);
xlabel('\omega (rad/s)'); ylabel('\angle H(j\omega) (deg)')

%% Print results
fprintf('Optimal sets found: \n\n')
fprintf('    Kp    tau    Wnm    Znm     Tl    Ti \n');
fprintf('A %-.4f %-.4f %-.4f %-.4f \n', Kp_A, tau_A, Wnm_A, Znm_A);
fprintf('B %-.4f %-.4f %-.4f %-.4f %-.4f \n', Kp_B, tau_B, Wnm_B, Znm_B,Tl_B);
fprintf('C %-.4f %-.4f %-.4f %-.4f %-.4f %-.4f \n', Kp_C, tau_C, Wnm_C, Znm_C, Tl_C, Ti_C);

save('done_part3.mat')