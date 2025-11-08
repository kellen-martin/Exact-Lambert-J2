function [v_1,v_2, n] = UV_Lambert_Bisect(r_1_vec, r_2_vec, mu, delta_t, tol, longway)
% UV_Lambert_Bisect solves Lambert's Problem in Universal Variables
% Using bisection method

    % Calculate geometry
    r_1 = norm(r_1_vec);
    r_2 = norm(r_2_vec);
    c_theta = (dot(r_1_vec, r_2_vec)/(r_1*r_2));
    if longway
        theta = (dot(r_1_vec, r_2_vec)/(r_1*r_2));
        theta = 2*pi - theta;
        c_theta = cos(theta);
    end
    A = sqrt(r_1*r_2*(1 + c_theta));

    % initial conditions for iteration
    psi = 0;
    c_2 = 1/2;
    c_3 = 1/6;

    psi_upper = 4*pi^2;
    psi_lower = -4*pi^2;

    n_max = 10000;

    % Bistection
    f_psi = 10000;
    n=0;
    while(abs(f_psi)>tol && n<=n_max)
        y = r_1 + r_2 + (A*(psi*c_3 - 1)/sqrt(c_2));
        chi = sqrt(y/c_2);

        f_psi = (chi^3*c_3 + A*sqrt(y))/sqrt(mu) - delta_t;

        if(f_psi<0)
            psi_lower = psi;
        else
            psi_upper = psi;
        end

        psi = (psi_upper + psi_lower)/2;
        c_2 = stumpff_c2(psi);
        c_3 = stumpff_c3(psi);

        n = n + 1;
    end
    if(n >= n_max)
        fprintf('did not converge')
    end
    % Calculate Lagrange Coefficients
    f = 1 - y/r_1;
    g = A*sqrt(y/mu);
    g_dot = 1 - y/r_2;

    % Calculate Velocities
    v_1 = (r_2_vec - f*r_1_vec)/g;
    v_2 = (g_dot*r_2_vec - r_1_vec)/g;

end