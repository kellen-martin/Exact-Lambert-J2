%% This Scipt is for testing the Generalized Labert's Problem 
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

% Time step
t_step = 1;        % [s]


%% Generate random orbit 
check = true;
n = 1;
N = 16;

while check
[r1, v1, r2, v2, delta_t, a_L] = lambert_conditions(mu);
% delta_t = 5.301332507081173e+03;
c = norm(r1 - r2);
semi = .5*(norm(r1) + norm(r2) + c);
a_min=semi/2;
theta = acos(dot(r1/norm(r1), r2/norm(r2)));
check = cross(r1, r2);

if(check(3)<0)
    theta = 2*pi-theta;
end



[a_test, e, p, inc, Omega1, ~, f1] = get_oe(r1, v1, mu);


[a2, ~, ~, i, Omega2, ~, f2] = get_oe(r2, v2, mu);
T = 2*pi*sqrt((a2^3)/mu);
ratio = delta_t/T;

[a, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);


a_error = (abs(a_test - a)/a_test);


% 
% error = norm(r2 - r2_sim)/norm(r2);
error_v2_rel = norm(v2_L - v2)/norm(v2);
error_v1_rel = norm(v1_L - v1)/norm(v1);


h1_hat = cross(r1, v1)/norm(cross(r1, v1));

gamma1 = acos(dot(r1, v1)/(norm(r1)*norm(v1)));
gamma2 = acos(dot(r2, v2)/(norm(r2)*norm(v2)));

% Propegate
[r2_check, ~] = pkepler(r1, v1_L, delta_t, 0, 0);
abs_err = norm(r2_check - r2);
rel_err = abs_err/norm(r2);
if rel_err > .02
    check = false;
end
n = n+1;
end

[a, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);
[at, et, pt, it, Omegat, omegat, ft] = get_oe(r1, v1_L, mu);

t_vec = linspace(0, delta_t, 100);
for i = 1:length(t_vec)
    [orb(:,i), ~] = pkepler(r1, v1, t_vec(i), 0, 0);
    [orb_sol(:,i), ~] = pkepler(r1, v1_L, t_vec(i), 0, 0);
end

figure
plot3(orb(1,:), orb(2,:), orb(3,:))
hold on
plot3(orb_sol(1,:), orb_sol(2,:), orb_sol(3,:))