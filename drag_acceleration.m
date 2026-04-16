function ap = drag_acceleration(zeta, v, bc)
H = 50;
rho0 = 1e-8;
zeta0 = 100;
rho = rho0*exp(-(zeta-zeta0)/H);
ap = (-0.5 ./ bc) .* rho .* v.^2;
end