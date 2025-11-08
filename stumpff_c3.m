function c3 = stumpff_c3(psi)
    if psi > 1e-6
        c3 = (sqrt(psi) - sin(sqrt(psi))) / (psi^(1.5));
    elseif psi < -1e-6
        c3 = (sinh(sqrt(-psi)) - sqrt(-psi)) / ((-psi)^(1.5));
    else
        c3 = 1/6;
    end
end