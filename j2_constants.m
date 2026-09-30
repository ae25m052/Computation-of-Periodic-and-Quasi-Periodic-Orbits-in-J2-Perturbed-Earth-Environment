function c = j2_constants()
% Earth constants, canonical units and integrator settings.

c.mu = 398600.4418;        % km^3/s^2
c.Re = 6378.1363;          % km
c.J2 = 1.08262668e-3;

c.DU = c.Re;
c.TU = sqrt(c.DU^3/c.mu);
c.VU = c.DU/c.TU;

c.opts = odeset('RelTol',1e-13,'AbsTol',1e-13);

end
