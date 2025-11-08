function [oes, t] = oeIntegrateJ2(oe_0, delta_t)
mu = 3.986*10^5;    % [km^3/s^2]

a0 = oe_0(1);
e0 = oe_0(2);
i0 = oe_0(3);
omega0 = oe_0(4);
Omega0 = oe_0(5);
f0 = oe_0(6);


h0 = sqrt(mu*a0*(1 - e0^2));

State0 = [h0; e0; f0; Omega0; i0; omega0];
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
tspan = [0 delta_t];

[t, State] = ode45(@(t, State) ode_func(State), tspan, State0, options);

h = State(:, 1);
e = State(:, 2);
f = State(:, 3);
Omega = State(:, 4);
i = State(:, 5);
omega = State(:, 6);

p = h.^2/mu;
a = p./(1 - e.^2);

oes = [a, e, i, omega, Omega, f];


end


function dStatedt = ode_func(State)

mu = 3.986*10^5;    % [km^3/s^2]
J_2 = 1.0826E-3;
Re = 6378;          % [km]

h = State(1);
e = State(2);
f = State(3);
Omega = State(4);
i = State(5);
omega = State(6);

% calculate radius
p = h^2/mu;
r = p/(1 + e*cos(f));
% From Curtis (pg 685)
u = omega + f;

dh = -(3/2)*J_2*mu*Re^2/(r^3)*sin(i)^2*sin(2*u);
de = (3/2)*J_2*mu*Re^2/(h*r^3)*(h^2/(mu*r)*sin(f)*(3*sin(i)^2*sin(u)^2 - 1) - sin(2*u)*sin(i)^2*((2 + e*cos(f))*cos(f) + e));
df = h/(r^2) + (3/2)*J_2*mu*Re^2/(e*h*r^3)*(h^2/(mu*r)*cos(f)*(3*sin(i)^2*sin(u)^2 - 1) + (2 + e*cos(f))*sin(2*u)*sin(i)^2*sin(f));
dOmega = -3*J_2*mu*Re^2/(h*r^3)*sin(u)^2*cos(i);
di = -(3/4)*J_2*mu*Re^2/(h*r^3)*sin(2*u)*sin(2*i);
domega = (3/2)*J_2*mu*Re^2/(e*h*r^3)*(h^2/(mu*r)*cos(f)*(1 - 3*sin(i)^2*sin(u)^2) - (2 + e*cos(f))*sin(2*u)*sin(i)^2*sin(f) + 2*e*cos(i)^2*sin(u)^2);

dStatedt = [dh; de; df; dOmega; di; domega];
end
