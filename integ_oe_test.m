%% Test of Classical Orbit Element Integration using ode45
clc
clear
close all

mu = 3.986*10^5;    % [km^3/s^2]
J_2 = 1.0826E-3;
alpha = 6378;          % [km]
N = 256;

[r1_lc, v1_lc, ~, ~, delta_t, ~] = lambert_conditions_lopez();

[a, e, p, i, Omega, omega, f] = get_oe(r1_lc, v1_lc, mu);
oes_0 = [a, e, i, omega, Omega, f];

delta_t = round(delta_t);
oes = oeIntegrateJ2(oes_0, delta_t);


for i = 1:length(oes)
    [r(i, :), ~] = get_rv(oes(i,:), mu);
end

[r1, v1] = get_rv(oes_0, mu);
[r2, v2] = get_rv(oes(end, :), mu);

[r2_L, v2_L, a_out] = lopezPropegate(r1, v1, delta_t);
[r2_V, ~] = pkepler(r1, v1, delta_t, 0, 0);

diff_L = norm(r2 - r2_L);
diff_V = norm(r2 - r2_V);

[a, v1_s, v2_s] = Lamabert_J2_1newt(r1_lc, r2, delta_t, mu, J_2, alpha, N);

%% Check with integ
[a, e, p, i, Omega, omega, f] = get_oe(r1, v1_s, mu);
oes_0s = [a, e, i, omega, Omega, f];

oes_s = oeIntegrateJ2(oes_0s, delta_t);

[r2_i_check, ~] = get_rv(oes_s(end, :), mu);
err = norm(r2_i_check - r2);

[r2_Ls, ~, ~] = lopezPropegate(r1, v1_s, delta_t);
err_L = norm(r2_Ls - r2);
