function [r, v] = coe2rv(a, e, Omega, i, omega, theta, mu)

p = a*(1-e^2);
r_pf = p/(1+e*cos(theta)) * [cos(theta); sin(theta); 0];
v_pf =  sqrt(mu/p) * [-sin(theta); e + cos(theta); 0];

R3_Omega = [cos(Omega) sin(Omega) 0;
    -sin(Omega) cos(Omega) 0;
    0 0 1];

R1_i = [1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];

R3_omega = [cos(omega) sin(omega) 0;
    -sin(omega) cos(omega) 0;
    0 0 1];

dcm = R3_Omega * R1_i * R3_omega;

r = dcm * r_pf;
v = dcm * v_pf;

end
