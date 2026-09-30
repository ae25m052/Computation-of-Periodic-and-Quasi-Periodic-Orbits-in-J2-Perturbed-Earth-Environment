function [r,v] = kep2rv(a,e,i,RAAN,argp,nu,mu)
% Classical elements to position and velocity. Angles in radians.

p    = a*(1-e^2);
rmag = p/(1+e*cos(nu));

r_pf = [rmag*cos(nu); rmag*sin(nu); 0];
v_pf = sqrt(mu/p)*[-sin(nu); e+cos(nu); 0];

R3_W = [ cos(-RAAN) sin(-RAAN) 0; -sin(-RAAN) cos(-RAAN) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)];
R3_w = [ cos(-argp) sin(-argp) 0; -sin(-argp) cos(-argp) 0; 0 0 1];

Q = R3_W*R1_i*R3_w;

r = Q*r_pf;
v = Q*v_pf;

end
