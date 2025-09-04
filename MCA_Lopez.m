%% This Scipt is for MCA testing against Vallado 
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
trials = 1;
a_error = zeros(1, trials);
abs_error = zeros(1, trials);
rel_error = zeros(1, trials);
%h = waitbar(0, 'Running');
% boxFolder = 'C:\Users\kmartin6\Box\';
% saveFile = fullfile(boxFolder, 'lopez_results.mat');
% cleanupObj = onCleanup(@() save(saveFile));

for i=1:trials
    % Random condition
    [r1, v1, r2, v2, delta_t, a_L] = lambert_conditions_lopez();
    [a_test, e, p, inc, Omega1, ~, f1] = get_oe(r1, v1, mu);
    
    % Solve
    [a, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);

    % Get Errors
    a_error(i) = (abs(a_test - a)/a_test);
    [r2_check, ~, ~] = lopezPropegate(r1, v1, delta_t);
    abs_error(i) = norm(r2_check - r2);
    rel_error(i) = abs_error(i)/norm(r2);
   % waitbar(i/trials, h)
end


%save(saveFile, 'abs_error', 'rel_error');
%close(h)
%% Post Process Data
abs_error_mean = mean(abs_error);
mu_a = std(abs_error);
rel_err_mean = mean(rel_error);
mu_r = std(rel_error);

%rel_error = rmoutliers(rel_error);
figure
histogram(rel_error, 100, 'Normalization','probability')
hold on
xline(rel_err_mean)
hold off
title('Relative Error Distribution')
xlabel('Relative Error')

%abs_error = rmoutliers(abs_error);
figure
histogram(abs_error, 100, 'Normalization','probability')
hold on
xline(abs_error_mean)
hold off
title('Absolute Error Distribution')
xlabel('Absolute Error')

