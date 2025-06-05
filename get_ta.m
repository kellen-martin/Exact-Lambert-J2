function [f1, f2] = get_ta(r1, r2, theta, p, e)
% Find the true anomaly (angle ambiguity)
f1_p = real(acos((p/norm(r1) - 1)/e));
f2_p = real(acos((p/norm(r2) - 1)/e));

% Test all combinations of f1, f2 to match theta
f1_pm = [f1_p, -f1_p, -f1_p, f1_p];
f2_pm = [f2_p, -f2_p, f2_p, -f2_p];
theta_check = mod(f2_pm - f1_pm, 2*pi);
diff = abs(theta - theta_check);
[~, index] = min(diff);

% Return the values of true anomaly
f1 = f1_pm(index);
f2 = f2_pm(index);
end