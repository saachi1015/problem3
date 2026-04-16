function ap = J2_acceleration(r)
Omega = 195 * pi/180; %rad
i = 54 * pi/180; %rad
w = 235 * pi/180; %rad
ta = 20 * pi/180; %rad
J2 = 0.000108263;

muE = 3.986e5;
rE = 6371;

function mat = R3(angle)
    mat = [cos(angle), sin(angle), 0;
        -sin(angle), cos(angle), 0;
        0, 0, 1];
end

function mat = R1(angle)
    mat = [1, 0, 0;
        0, cos(angle), sin(angle);
        0, -sin(angle), cos(angle)];
end

dcm = R3(Omega)*R1(i)*R3(w);
ur = [cos(ta); sin(ta); 0];

u_ECI = dcm * ur;

ap = (3/2)*J2* (muE*(rE)^2)./(r.^4) * ((5*(r(3)^2)/(norm(r)^2)) -1)*u_ECI.' - (2* (r(3)/ norm(r))*[0,0,1]);

end