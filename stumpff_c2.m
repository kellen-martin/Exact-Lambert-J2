function c2 = stumpff_c2(psi)
    if psi > 1e-6
        c2 = (1 - cos(sqrt(psi))) / psi;
    elseif psi < -1e-6
        c2 = (1 - cosh(sqrt(-psi))) / psi;
    else
        c2 = 0.5;
    end
end