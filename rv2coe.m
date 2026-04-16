function [a,e,i,Omega,omega,theta_ta,e_vec,h,n] = rv2coe(r0,v0)

    % constants
    muE = 3.986e5;  % km^3/s^2

    % basis vectors
    xhat = [1;0;0];
    zhat = [0;0;1];

    % make sure inputs are columns
    r0 = r0(:);
    v0 = v0(:);

    % magnitudes
    r = norm(r0);
    v = norm(v0);

    % angular momentum
    h = cross(r0, v0);
    hmag = norm(h);

    % node vector
    n = cross(zhat, h);
    nmag = norm(n);

    % eccentricity vector
    e_vec = (cross(v0, h) / muE) - (r0 / r);
    e = norm(e_vec);

    % specific orbital energy
    eps = v^2/2 - muE/r;

    % semi-major axis
    if abs(e - 1) > 1e-12
        a = -muE / (2*eps);
    else
        a = Inf; % parabolic case
    end

    % inclination
    i = acos(h(3)/hmag);

    % RAAN
    if nmag > 1e-12
        Omega = acos(n(1)/nmag);
        if n(2) < 0
            Omega = 2*pi - Omega;
        end
    else
        Omega = 0; % equatorial orbit
    end

    % argument of periapsis
    if nmag > 1e-12 && e > 1e-12
        omega = acos(dot(n, e_vec)/(nmag*e));
        if e_vec(3) < 0
            omega = 2*pi - omega;
        end
    else
        omega = 0; % circular and/or equatorial special case
    end

    % true anomaly
    if e > 1e-12
        theta_ta = acos(dot(e_vec, r0)/(e*r));
        if dot(r0, v0) < 0
            theta_ta = 2*pi - theta_ta;
        end
    else
        theta_ta = 0; % circular special case
    end
end