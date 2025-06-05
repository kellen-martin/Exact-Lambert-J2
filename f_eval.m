function val = f_eval(theta, mu, J2, alpha, c, a, semi, r1, r2, delta_t)
    % find the angles
    [psi, cphi] = lagrange_angles(semi, c, a, theta);
    e = get_eccentricity(r1, r2, theta, psi, a);
    [i,~,~] = newton_angles(r1, r2, a, e, J2, mu, alpha, delta_t);
    val = (2*(psi - sin(psi)*cphi) + (J2*((alpha/(2*a))^2)*(3*sin(i)^2 - 2))/(((1 - e^2)^3))*(4*(e^2 + 2)*psi - 16*sin(psi)*cphi + 8*sin(2*psi)*cphi^2 - 4*e^2*sin(2*psi)) - sqrt(mu/(a^3))*delta_t)^-1;
end