function e = get_eccentricity(r1_vec, r2_vec, theta, psi, a)
    r1 = norm(r1_vec);
    r2 = norm(r2_vec);
    
    % Calculate Orbit Parameter
    p = norm(r1)*norm(r2)*sin(theta/2)^2/(a*sin(psi)^2);

    % Claculate eccentricity
    e = sqrt(1 - p/a);
end
