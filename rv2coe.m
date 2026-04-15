function [a, e, omega, i, w, theta_ta] = rv2coe(r0,v0)
% define x and z directions in ECI
x = [1;0;0];
z = [0;0;1];

muE = 3.986e5; % constant
r = norm(r0);
h = cross(r0, v0); % angular momentum
e_vec = (cross(v0, h) / muE) - (r0 / r); % eccentricity vector
e = norm(e_vec); % eccentricity value
a = dot(h,h) / (muE*(1-e^2)); % semi-major axis
i = mod(acos(dot(h,z)/norm(h)), 2*pi); % angle of inclination
n = cross(z,h).';% line of nodes
omega = mod(acos(dot(n,x)/(norm(x)*norm(n))), 2*pi); % right ascension of the
% ascending node
w = mod(acos(dot(n,e_vec)/(norm(n)*e)),2*pi); % argument of perapsis
theta_ta = mod(acos(dot(e_vec,r0)/(e*r)),2*pi); % true anomaly

end