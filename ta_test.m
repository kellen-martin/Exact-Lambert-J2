%% Test of finding true anomalies
clc
clear
%% System Variables
% Gravitational parameter
mu = 3.986*10^5;    % [km^3/s^2]

% Value of J_2
J_2 = 1.0826E-3;

% Mean Equitorial Radius
alpha = 6378;       % [km]

% Time step
t_step = 1;        % [s]


%% Generate random orbit 

n = 0;
N = 128;

[r1, v1, r2, v2, delta_t, a_L] = lambert_conditions(mu);


[a_test, e, p, inc, Omega1, ~, f1] = get_oe(r1, v1, mu);
[a2, ~, ~, i, Omega2, ~, f2] = get_oe(r2, v2, mu);

theta = acos(dot(r1/norm(r1), r2/norm(r2)));

check = cross(r1, r2);

if(check(3)<0)
    theta = 2*pi-theta;
end

[f1_t, f2_t] = get_ta(r1, r2, theta, p, e);