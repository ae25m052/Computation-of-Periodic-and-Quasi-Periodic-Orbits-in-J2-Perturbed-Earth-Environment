function ref = eq_reference_apsides(h, e, mu, J2, Re)
%EQ_REFERENCE_APSIDES  Independent (non-shooting) check for eq_dc_pseudo.m.
%
%   ref = EQ_REFERENCE_APSIDES(h, e, mu, J2, Re)
%
%   Uses only the two conserved quantities, energy E and angular momentum h
%   (Wang et al. 2014, eqs. (4)-(8)). No orbit is integrated.
%
%   1. Apsides. At r_min and r_max the radial speed is zero, so both have
%      the same effective potential
%          Veff(r) = h^2/(2 r^2) - mu/r - J2 mu Re^2/(2 r^3)
%      Solve Veff(r_min) = Veff(r_min/rho), with rho = (1-e)/(1+e).
%   2. Energy E = Veff(r_min).
%   3. Third root r* of eq. (9):  E r^3 + mu r^2 - h^2 r/2 + J2 mu Re^2/2 = 0.
%      The product of the three roots is -(J2 mu Re^2/2)/E.
%   4. Radial period and angle per radial period, by quadrature of eq. (8):
%          T_r    = 2 * integral  dr / rdot
%          dtheta = 2 * integral (h/r^2) dr / rdot
%      Substituting r = (r_min+r_max)/2 - (r_max-r_min)/2*cos(phi) removes
%      the 1/0 at both ends, leaving the smooth integrand
%          dt/dphi = sqrt( r^3 / (2 (-E) (r - r*)) ),   phi from 0 to pi.
%
%   This is a different route from the paper's elliptic-integral formula
%   (eq. 16). A separate Python check found the two agree to ~1e-13 for all
%   five Table 1 orbits. That formula needs the elliptic integral of the
%   third kind, which base MATLAB (no toolboxes) does not have, so it is
%   not used here.

rho  = (1 - e)/(1 + e);
Veff = @(r) h^2./(2*r.^2) - mu./r - J2*mu*Re^2./(2*r.^3);

rc   = (h^2 + sqrt(h^4 - 6*J2*mu^2*Re^2))/(2*mu);    % bottom of the well
fz   = optimset('TolX', 1e-15);
rmin = fzero(@(r) Veff(r) - Veff(r/rho), [rc*rho, rc], fz);
rmax = rmin/rho;
E    = Veff(rmin);
rs   = -(J2*mu*Re^2/2)/(E*rmin*rmax);

rr     = @(p) (rmin + rmax)/2 - (rmax - rmin)/2*cos(p);
dtdphi = @(p) sqrt(rr(p).^3 ./ (2*(-E)*(rr(p) - rs)));
Tr     = 2*integral(dtdphi, 0, pi, 'AbsTol', 1e-14, 'RelTol', 1e-13);
dth    = 2*integral(@(p) h./rr(p).^2 .* dtdphi(p), 0, pi, ...
                    'AbsTol', 1e-14, 'RelTol', 1e-13);

ref.rmin   = rmin;
ref.rmax   = rmax;
ref.E      = E;
ref.rstar  = rs;
ref.Tr     = Tr;
ref.dtheta = dth;
end
