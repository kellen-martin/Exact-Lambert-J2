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
% N = [4 6 8 16 32 64];
% a_error_mean = zeros(1, length(N));
% v1_error_mean = zeros(1, length(N));
% v2_error_mean = zeros(1, length(N));
% solve_time_mean = zeros(1, length(N));
% 
% for j=1:length(N)
% error_v1_rel = zeros(1,1000);
% error_v2_rel = zeros(1,1000);
% error_a_rel = zeros(1,1000);
% solve_time = zeros(1,1000);
% delta_t = zeros(1,1000);

    for i=1:100
    % Create Random Initial Conditions
    
    [r1, v1, r2, v2, delta_t(i)] = lambert_conditions(mu);
    [a_v, e_v, p_v, i_v, Omega_1_v, ~, ~] = get_oe(r1, v1, mu);
    [i_1, Omega_11, Omega_21] = newton_angles(r1, r2, a_v, e_v, J_2, mu, alpha, delta_t(i));
    
    T = 2*pi*sqrt((a_v^3)/mu);
    ratio = delta_t/T;
    
    % Lambert Solve
    tic
    [a, v1_mag, v2_mag] = Lamabert_J2(r1, r2, delta_t(i), mu, J_2, alpha, N(j));
    solve_time(i) = toc;
    
    % Errors
    v1_mag_real = norm(v1);
    v2_mag_real = norm(v2);
    error_v1_rel(i) = abs(v1_mag - v1_mag_real)/v1_mag_real;
    error_v2_rel(i) = abs(v2_mag - v2_mag_real)/v2_mag_real;
    error_a_rel(i) = abs(a - a_v)/abs(a_v);
    end

% Mean Error
a_error_mean(j) = mean(error_a_rel);
v1_error_mean(j) = mean(error_v1_rel);
v2_error_mean(j) = mean(error_v2_rel);
solve_time_mean(j) = mean(solve_time);
% end
%% Plots
run_time_analysis = table(N', solve_time_mean', a_error_mean');


% Error Plots
figure
histogram(error_a_rel, 100)
hold on
xline(a_error_mean, '--r', 'LineWidth',1.5)
xlim ([0, max(error_a_rel)]);
title('Relative Error in Semi-major Axis')
xlabel('Relative Error')
ylabel('Number of Trials')
legend('Relative Errors', 'Mean')

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