%% Test of Classical Orbit Element Integration using ode45
mu = 3.986*10^5;    % [km^3/s^2]
J_2 = 1.0826E-3;
alpha = 6378;          % [km]
N = 256;

[oes_0, delta_t] = random_orbit(mu);
oes = oeIntegrateJ2(oes_0, delta_t);

[r1, v1] = get_rv(oes(1, :), mu);
[r2, v2] = get_rv(oes(end, :), mu);

[r2_L, v2_L, a_out] = lopezPropegate(r1, v1, delta_t);
[r2_V, ~] = pkepler(r1, v1, delta_t, 0, 0);

diff_L = norm(r2 - r2_L);
diff_V = norm(r2 - r2_V);

[a, v1_s, v2_s] = Lamabert_J2_1newt(r1, r2, delta_t, mu, J_2, alpha, N);

%% Check with integ
[a, e, p, i, Omega, omega, f] = get_oe(r1, v1_s, mu);
oes_0s = [a, e, i, omega, Omega, f];

oes_s = oeIntegrateJ2(oes_0s, delta_t);

[r2_i_check, ~] = get_rv(oes_s(end, :), mu);
