function e = get_eccentricity(r1_vec, r2_vec, theta, psi, a)
    % Calculate Orbit Parameter
    p = get_parameter(r1_vec, r2_vec, theta, psi, a);

    % Claculate eccentricity
    e = sqrt(1 - p/a);
end
