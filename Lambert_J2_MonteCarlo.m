%% Monte Carlo Analysis
clc
clear

%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

%% Monte Carlo
N = 500;
trials = 100;
error_v1_rel = zeros(1,trials);
error_v2_rel = zeros(1,trials);
error_a_rel = zeros(1,trials);
delta_t = zeros(1,trials);
oes = zeros(6,trials);
h = waitbar(0, 'Running');
    for i=1:trials
    % Create Random Initial Conditions
    
    [r1, v1, r2, v2, delta_t(i), a_L] = lambert_conditions_lopez();
    [a2(i), e2(i), p2(i), i2(i), Omega2(i), omega2(i), f2(i)] = get_oe(r2, v2, mu);
    
    % Lambert Solve
    [a, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t(i), mu, J_2, alpha, N);

    % Errors
    error_a_rel(i) = abs(a_L - a)/a_L;
    error_v1_rel(i) = norm(v1_L - v1)/norm(v1);
    error_v2_rel(i) = norm(v2_L - v2)/norm(v2);
    waitbar(i/trials, h)
    end

% Mean Error
a_error_mean = mean(error_a_rel);
v1_error_mean = mean(error_v1_rel);
v2_error_mean = mean(error_v2_rel);

% Standard Deviation of the Error
v1_error_std = std(error_v1_rel);
v2_error_std = std(error_v2_rel);

%% Plots

% Error Plots
% figure
% histogram(error_a_rel, 100)
% hold on
% xline(a_error_mean, '--r', 'LineWidth',1.5)
% xlim ([0, max(error_a_rel)]);
% title('Relative Error in Semi-major Axis')
% xlabel('Relative Error')
% ylabel('Number of Trials')
% legend('Relative Errors', 'Mean')

figure
histogram(error_v1_rel, 100)
hold on
xline(v1_error_mean, '--r', 'LineWidth',1.5)
xlim ([0, max(error_v1_rel)]);
title('Relative Error in Initial Velocity')
xlabel('Relative Error')
ylabel('Number of Trials')
legend('Relative Errors', 'Mean')

figure
histogram(error_v2_rel, 100)
hold on
xline(v2_error_mean, '--r', 'LineWidth',1.5)
xlim ([0, max(error_v2_rel)]);
title('Relative Error in Final Velocity')
xlabel('Relative Error')
ylabel('Number of Trials')
legend('Relative Errors', 'Mean')