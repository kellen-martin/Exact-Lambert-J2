function [r1, v1, r2, v2, delta_t, oes] = lambert_conditions(mu)

% Random Initial Contitiobns
[oes, delta_t] = random_orbit(mu);
[r1, v1] = get_rv(oes, mu);

% Use Vallado perturbed kepler solver to find final condition
[r2,v2] = pkepler(r1, v1, delta_t, 0, 0);
end