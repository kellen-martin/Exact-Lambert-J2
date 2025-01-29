function [i, Omega_1, Omega_2] = newton_angles(r1_vec, r2_vec, a, e, J2, mu, alpha, delta_t)
    Beta = 3/2*((J2*sqrt(mu)*alpha^2)/(a^(7/2)*(1 - e^2)^2));

    F = @(x) [sin(x(1))*sin(x(2))*r1_vec(1) - sin(x(1))*cos(x(2))*r1_vec(2) + cos(x(1))*r1_vec(3);
              sin(x(1))*sin(x(3))*r2_vec(1) - sin(x(1))*cos(x(3))*r2_vec(2) + cos(x(1))*r2_vec(3);
              x(3) + Beta*cos(x(1))*delta_t - x(2)];

    % Make better initial guess for inclinaction
    options = optimset('Display','off');
    options.MaxFunEvals=100000;
    options.MaxIter = 100000;
    x = fsolve(F, [pi/4, 0, 0], options);

    r_cross_r = cross(r1_vec, r2_vec)/norm(cross(r1_vec, r2_vec));

    i = x(1);
    Omega_1 = x(2);
    Omega_2 = x(3);

    h = [sin(x(1))*sin(x(2)), -sin(x(1))*cos(x(2)), cos(x(1))];

    if (r_cross_r(3)*h(3) < 0)
        i = i + pi;
    end
end