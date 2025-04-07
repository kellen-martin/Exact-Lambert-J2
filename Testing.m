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
n = 0;
N = 128;

while(check)
[r1, v1, r2, v2, delta_t, a_L] = lambert_conditions(mu);
% r1 = [-2.010663343027483e+03;-7.711567092820798e+03;-2.826462291832482e+03];
% r2 = [1.411398563499769e+04;-1.006592216200454e+03;-7.083400698437274e+03];
% v1 = [4.807826587768576;-3.977024328378902;-4.195440314892678];
% v2 = [0.006443467663318;3.192992366201820;1.561820415609631];
% delta_t = 5.301332507081173e+03;
c = norm(r1 - r2);
semi = .5*(norm(r1) + norm(r2) + c);
a_min=semi/2;
theta = acos(dot(r1/norm(r1), r2/norm(r2)));
% check = cross(r1, r2);
% 
% if(check(3)<0)
%     theta = 2*pi-theta;



[a_test, e, p, inc, Omega1, ~, f1] = get_oe(r1, v1, mu);
[a2, ~, ~, i, Omega2, ~, f2] = get_oe(r2, v2, mu);
[v1_test, v2_test] = velocity_solve(r1, r2, norm(v1), norm(v2), Omega1, Omega2, inc, p, mu);

[a, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);

a_error = (abs(a_test - a)/a_test);

[r2_sim,v2_sim] = pkepler(r1, v1_L, delta_t, 0, 0);
% 
% error = norm(r2 - r2_sim)/norm(r2);
error_v2_rel = norm(v2_L - v2)/norm(v2);
error_v1_rel = norm(v1_L - v1)/norm(v1);
error_v2_test = norm(v2_test - v2)/norm(v2);

check = error_v2_rel < .2;
n = n+1;
end

gama_2 = asin(norm(cross(r2,v2))/(norm(r2)*norm(v2)));
r2_hat = r2/norm(r2);
h = cross(r2, v2);
h2_hat = h/norm(h);
theta2_hat = cross(h2_hat, r2_hat);
v2_aaa = norm(v2)*(-cos(gama_2)*r2_hat + sin(gama_2)*theta2_hat);

[a, v1_L, v2_L] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);

val1 = f_check(theta, mu, J_2, alpha, c, a2, semi, r1, r2, delta_t);
val2 = f_check(theta, mu, J_2, alpha, c, a, semi, r1, r2, delta_t);

