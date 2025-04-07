function [oes, delta_t] = random_orbit(mu)
% Creates a random orbit, bounds are given in the function, may change this
% later

a_min = 6400;
a_max = 11000;
a = a_min + (a_max - a_min)*rand;

e_max = .8;
e = e_max*rand;

i_max = pi/2;
i = i_max*rand;

omega = pi*rand;

Omega = pi*rand;

f = 2*pi*rand;

oes = [a, e, i, omega, Omega, f];

% Currently bounded by half the orbit period
T_max = pi*sqrt((a^3)/mu);
delta_t = T_max*rand;

end