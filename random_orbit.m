function [oes, delta_t] = random_orbit(mu)
% Creates a random orbit, bounds are given in the function, may change this
% later
rp_min = 6400;
a_max = 11000;

e_max = .5;
e = e_max*rand;

a_min = rp_min/(1-e);
a = a_min + (a_max - a_min)*rand;

i_max = pi/2;
i = i_max*rand;

omega = pi*rand;

Omega = pi*rand;

f = 2*pi*rand;
oes = [a, e, i, omega, Omega, f];

% Currently bounded by half the orbit period
t_min = 100;
T_max = 2*pi*sqrt((a^3)/mu) -100;
delta_t = t_min + (T_max - t_min)*rand;

end