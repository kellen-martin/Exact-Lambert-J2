function dStatedt = ode_fun_J2(State, mu)
r = State(1:3);
v = State(4:6);
rdot = v;

% J2 Stuff 
J_2 = 1.0826E-3;
R = 6378; % [km]
vec = [(r(1)^2 + r(2)^2 - 4*r(3)^2)*r(1); 
     (r(1)^2 + r(2)^2 - 4*r(3)^2)*r(2);
     (3*(r(1)^2 + r(2)^2) - 2*r(3)^2)*r(3)];
u_J2 = (-(3*mu/(2*norm(r)^7))*R^2*J_2^2)*vec;


vdot = -mu/(norm(r)^3).*r + u_J2;

dStatedt = [rdot; vdot];
end 