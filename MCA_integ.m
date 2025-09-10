%% This Scipt is for MCA testing against Integration of the OEs 
clc;
clear
close all

%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

N = 256;

%% MCA
trials = 500;
a_error = zeros(1, trials);
abs_error = zeros(1, trials);
rel_error = zeros(1, trials);
h = waitbar(0, 'Running');
for i=1:trials
    % Random condition
    [r1, v1, r2, v2, delta_t, a_L] = lambert_conditions(mu);
    % [a, e, p, inc, Omega, omega, f] = get_oe(r1, v1, mu);
    
    % Solve
    [~, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);
    [a, e, p, inc, Omega, omega, f] = get_oe(r1, v1_L, mu);
    oes_0 = [a, e, inc, omega, Omega, f];

    % Get Errors
    oes = oeIntegrateJ2(oes_0, delta_t);
    [r2_check, ~] = get_rv(oes(end, :), mu);

    abs_error(i) = norm(r2_check - r2);
    rel_error(i) = abs_error(i)/norm(r2);
    waitbar(i/trials, h)
end
close(h)
%% Post Process Data
% abs_error_mean = mean(abs_error);
% mu_a = std(abs_error);
% rel_err_mean = mean(rel_error);
% mu_r = std(rel_error);
% 
% rel_error = rmoutliers(rel_error);
% figure
% histogram(rel_error, 100, 'Normalization','probability')
% hold on
% xline(rel_err_mean)
% hold off
% title('Relative Error Distribution')
% xlabel('Relative Error')
% 
% abs_error = rmoutliers(abs_error);
% figure
% histogram(abs_error, 100, 'Normalization','probability')
% hold on
% xline(abs_error_mean)
% hold off
% title('Absolute Error Distribution')
% xlabel('Absolute Error')

