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
trials = 10000;
a_error = zeros(1, trials);
abs_error = zeros(1, trials);
rel_error = zeros(1, trials);
abs_errorv = zeros(1, trials);
rel_errorv = zeros(1, trials);
tof = zeros(1, trials);
theta = zeros(1, trials);
ecc = zeros(1,trials);
a_min_close = zeros(1, trials);

h = waitbar(0, 'Running');
for i=1:trials
    % Random condition
    isBad = true;
    while(isBad)
        [r1, v1, r2, v2, delta_t, a_L] = lambert_conditions(mu);
        tof(i) = delta_t;
        c = norm(r1 - r2);
        semi = .5*(norm(r1) + norm(r2) + c);
        a_min=semi/2;
        theta(i) = acos(dot(r1/norm(r1), r2/norm(r2)));
        check = cross(r1, r2);
        
        longway = false;
        if(check(3)<0)
            theta(i) = 2*pi-theta(i);
            longway = false;
        end

        if (2.8 < theta(i) && theta(i)< 3.2) || (theta(i) > 5.8)
        else
            isBad = false;
        end
    end
    % Solve
    [~, v1_L, v2_L] = Lambert_J2_1newt_2(r1, r2, delta_t, mu, J_2, alpha, N);
    tol = 1e-9;
    [v_1v,v_2v, nv] = UV_Lambert_Bisect(r1, r2, mu, delta_t, tol, longway);
    [a, e, p, inc, Omega, omega, f] = get_oe(r1, v1_L, mu);
    [av, ev, pv, incv, Omegav, omegav, fv] = get_oe(r1, v_1v, mu);
    oes_0 = [a, e, inc, omega, Omega, f];
    oes_0v = [av, ev, incv, omegav, Omegav, fv];
    
    % Get Errors
    oes = oeIntegrateJ2(oes_0, delta_t);
    oesv = oeIntegrateJ2(oes_0v, delta_t);
    [r2_check, ~] = get_rv(oes(end, :), mu);
    [r2_checkv, ~] = get_rv(oesv(end, :), mu);
    % [r2_check, ~] = pkepler(r1, v1_L, delta_t, 0, 0);
    % [r2_checkv, ~] = pkepler(r1, v_1v, delta_t, 0, 0);


    ecc(i) = e;
    a_min_close(i) = (a - a_min)/a;
    abs_error(i) = norm(r2_check - r2);
    rel_error(i) = abs_error(i)/norm(r2);

    abs_errorv(i) = norm(r2_checkv - r2);
    rel_errorv(i) = abs_error(i)/norm(r2);

    waitbar(i/trials, h)
end
close(h)
%% Post Process Data
abs_error_mean = mean(abs_error);
abs_error_meanv = mean(abs_errorv);
mu_a = std(abs_error);
rel_err_mean = mean(rel_error);
mu_r = std(rel_error);

figure
histogram(rel_error, 100, 'Normalization','probability')
hold on
xline(rel_err_mean)
hold off
title('Relative Error Distribution')
xlabel('Relative Error')

figure
histogram(abs_error, 100, 'Normalization','probability')
hold on
xline(abs_error_mean)
hold off
title('Absolute Error Distribution')
xlabel('Absolute Error')

figure
scatter(theta, abs_error, '*')
title('Transfer Angle vs Error')
xlabel('Transfer Angle [rad]')
ylabel('Absolute Error [km]')

figure
scatter(tof, abs_error, '*')
title('Time-of-Flight vs Error')
xlabel('Time-of-Flight [sec]')
ylabel('Absolute Error [km]')

figure
scatter(a_min_close, abs_error, '*')
title('Error near a_{min}')
xlabel('Releative Difference between Solution and a_{min}')
ylabel('Absolute Error [km]')

figure
scatter(ecc, abs_error, '*')
title('Eccentricity vs Error')
xlabel('Eccentricity')
ylabel('Absolute Error [km]')

%% Remove near a-min
idx = find(a_min_close > .02);
err_rm = abs_error(idx);
err_rm_mean = mean(err_rm);
tof_rm = tof(idx);
theta_rm = theta(idx);

figure
histogram(err_rm, 100, 'Normalization','probability')
hold on
xline(err_rm_mean)
hold off
title('Absolute Error Distribution')
xlabel('Absolute Error')

figure
scatter(tof_rm, err_rm, '*')
title('Time-of-Flight vs Error')
xlabel('Time-of-Flight [sec]')
ylabel('Absolute Error [km]')

figure
scatter(theta_rm, err_rm, '*')
title('Transfer Angle vs Error')
xlabel('Transfer Angle [rad]')
ylabel('Absolute Error [km]')

figure
scatter(a_min_close(idx), err_rm, '*')
title('Error near a_{min}')
xlabel('Releative Difference between Solution and a_{min}')
ylabel('Absolute Error [km]')

%% Better Looking
err_hist = rmoutliers(err_rm);
figure
histogram(err_hist, 100, 'Normalization','probability')
hold on
xline(err_rm_mean)
hold off
title('Absolute Error Distribution')
xlabel('Absolute Error')