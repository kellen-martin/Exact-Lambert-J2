function [r, v] = integrate_2bp_j2(r0, v0, mu, time)
State0 = [r0; v0];

% Givens 
R_0 = 6378; % [km]
J_2 = 1.0826E-3; 

% ODE Call 
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
tspan = [0 time];
[t,State] = ode45(@(t, State) ode_fun_J2(State, mu), tspan, State0, options);

% Decompose State Vector
r = State(:, 1:3);
v = State(:, 4:6);

end