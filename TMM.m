clear;clc;

% --- 常量 ---
c0 = 299792458; % 真空中的光速 (m/s)
mu0 = 4 * pi * 1e-7; % 真空磁导率 (H/m)
eps0 = 1 / (mu0 * c0^2); % 真空介电常数 (F/m)

% --- 材料参数 ---
eps_r1 = 1; % 空气的相对介电常数 (epsilon_0)
eps_r2 = 4; % 介质的相对介电常数 (epsilon_1)
mu_r1 = 1; % 空气的相对磁导率
mu_r2 = 1; % 介质的相对磁导率

% --- 几何参数 ---
d = 10e-3; % 平板的厚度 (米), 10mm

% --- 频率范围 ---
freq_min = 1e9; % 1 GHz
freq_max = 20e9; % 20 GHz
num_points = 500;
f = linspace(freq_min, freq_max, num_points);
omega = 2 * pi * f;

% --- 本征阻抗 ---
eta1 = sqrt(mu_r1 * mu0 / (eps_r1 * eps0)); % 介质1（空气）的阻抗
eta2 = sqrt(mu_r2 * mu0 / (eps_r2 * eps0)); % 介质2（介质）的阻抗

% --- 波数 ---
k1 = omega * sqrt(mu_r1 * mu0 * eps_r1 * eps0); % 介质1（空气）中的波数
k2 = omega * sqrt(mu_r2 * mu0 * eps_r2 * eps0); % 介质2（介质）中的波数

% --- 单界面计算 ---
r_single = (eta2 - eta1) / (eta2 + eta1);
t_single = 2 * eta2 / (eta2 + eta1);

% --- 有限厚度平板的TMM计算 ---
r_slab = zeros(1, num_points);
t_slab = zeros(1, num_points);

for i = 1:num_points
    % 平板的传播矩阵 (P)，厚度为 d
    P = [exp(1i * k2(i) * d), 0; 0, exp(-1i * k2(i) * d)];

    % 介质1到介质2的界面矩阵 (I12)
    I12 = 0.5 * [1 + eta2/eta1, 1 - eta2/eta1; 1 - eta2/eta1, 1 + eta2/eta1];

    % 介质2到介质1的界面矩阵 (I21)
    I21 = 0.5 * [1 + eta1/eta2, 1 - eta1/eta2; 1 - eta1/eta2, 1 + eta1/eta2];

    % 总传输矩阵 (T)
    T = I12 * P * I21;

    % 从T矩阵中提取反射和透射系数
    r_slab(i) = - T(2, 1) / T(2, 2);
    t_slab(i) = T(1, 1) - T(1, 2)*T(2, 1)/ T(2, 2);
end

% --- 功率系数（复数的平方绝对值）---
R_single = abs(r_single)^2;
R_slab = abs(r_slab).^2;

T_single = 2 * abs(t_single)^2;
T_slab = abs(t_slab).^2;

% --- 绘图 ---
figure;
plot(f / 1e9, R_single * ones(size(f)), 'r--', 'LineWidth', 2);
hold on;
plot(f / 1e9, T_single * ones(size(f)), 'g--', 'LineWidth', 2);
plot(f / 1e9, R_slab, 'b-', 'LineWidth', 2);
plot(f / 1e9, T_slab, 'm-', 'LineWidth', 2);
grid on;
title('Reflection and transmission power coefficients');
xlabel('frequency (GHz)');
ylabel('power coefficients');
ylim([0, 1]);
legend('single-boundary reflection(|r|^2)', 'single-boundary transmission((n2/n1)*|t|^2)', 'slab reflection(|r|^2)', 'slab transmision(|t|^2)', 'Location', 'best');

figure;
plot(f / 1e9, abs(r_single) * ones(size(f)), 'r--', 'LineWidth', 2);
hold on;
plot(f / 1e9, abs(t_single) * ones(size(f)), 'g--', 'LineWidth', 2);
plot(f / 1e9, abs(r_slab), 'b-', 'LineWidth', 2);
plot(f / 1e9, abs(t_slab), 'm-', 'LineWidth', 2);
grid on;
title('Reflection and transmission coefficients');
xlabel('frequency (GHz)');
ylabel('coefficients');
ylim([0, 1]);
legend('single-boundary reflection(|r|)', 'single-boundary transmission(|t|)', 'slab reflection(|r|)', 'slab transmision(|t|)', 'Location', 'best');