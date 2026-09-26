function J = dxdoe(oe, mu)
%DXDOE  Analytical Jacobian d(r,v)/d(a,e,i,RAAN,argp,M).  6x6, closed form.
%
%   J = dxdoe(oe, mu)      oe = [a; e; i; RAAN; argp; M]   (M = mean anomaly)
%
%   Every entry is a hand-derived closed-form expression.  Nothing here is
%   numerically differenced.  Verified against central finite differences to
%   ~9e-10, which is the floor of that test at h = 1e-7.
%
%   Construction: in the perifocal frame, using eccentric anomaly E,
%
%       r_pf = [ a(cos E - e) ;  a*b*sin E ; 0 ],      b = sqrt(1-e^2)
%       v_pf = (eta/r)[ -sin E ; b*cos E ; 0 ],        eta = sqrt(mu*a)
%       r    = a(1 - e cos E),      M = E - e sin E
%
%   and r = R(RAAN,i,argp) r_pf.  The partials with respect to a, e and M act
%   on the perifocal vectors (through dE/de = a sinE/r and dE/dM = a/r); the
%   partials with respect to i, RAAN and argp act on the rotation matrix.
%
%   This is the object that makes an ANALYTICAL Keplerian STM possible --
%   see stm_kepler_analytic.m.

a = oe(1); e = oe(2); i = oe(3); Om = oe(4); w = oe(5); M = oe(6);

E = kepler_E(M, e);
s = sin(E); c = cos(E);
b = sqrt(1 - e^2);
r = a*(1 - e*c);
eta = sqrt(mu*a);

R = rot_pqw(Om, i, w);

p = [ a*(c - e) ; a*b*s ; 0 ];          % perifocal position
q = [ -eta*s/r  ; eta*b*c/r ; 0 ];      % perifocal velocity

E_e = a*s/r;        % dE/de at fixed M
E_M = a/r;          % dE/dM

% ---- d/da  (E unchanged at fixed e, M) --------------------------
p_a = [ c - e ; b*s ; 0 ];
q_a = -q/(2*a);                          % since q ~ a^{-1/2}

% ---- d/de  (at fixed M) -----------------------------------------
r_e = a*(-c + a*e*s^2/r);
p_e = [ -a*(1 + a*s^2/r) ; ...
         a*(-e*s/b + b*c*a*s/r) ; 0 ];
q_e = [ -eta*s*(a*c - r_e)/r^2 ; ...
         eta*((-e*c/b - a*b*s^2/r)*r - b*c*r_e)/r^2 ; 0 ];

% ---- d/dM --------------------------------------------------------
p_M = [ -a^2*s/r ; a^2*b*c/r ; 0 ];
q_M = [ -eta*(a*c - a^2*e*s^2/r)/r^2 ; ...
         eta*b*(-a*s - a^2*e*s*c/r)/r^2 ; 0 ];

J = zeros(6,6);
J(1:3,1) = R*p_a;              J(4:6,1) = R*q_a;
J(1:3,2) = R*p_e;              J(4:6,2) = R*q_e;
J(1:3,3) = drot_di(Om,i,w)*p;  J(4:6,3) = drot_di(Om,i,w)*q;
J(1:3,4) = drot_dOm(Om,i,w)*p; J(4:6,4) = drot_dOm(Om,i,w)*q;
J(1:3,5) = drot_dw(Om,i,w)*p;  J(4:6,5) = drot_dw(Om,i,w)*q;
J(1:3,6) = R*p_M;              J(4:6,6) = R*q_M;
end

% -----------------------------------------------------------------
function E = kepler_E(M, e)
E = M;  if e > 0.8, E = pi; end
for k = 1:100
    d = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - d;
    if abs(d) < 1e-15, break; end
end
end

function R = rot_pqw(Om, i, w)
cO=cos(Om); sO=sin(Om); ci=cos(i); si=sin(i); cw=cos(w); sw=sin(w);
R = [ cO*cw-sO*sw*ci , -cO*sw-sO*cw*ci ,  sO*si ; ...
      sO*cw+cO*sw*ci , -sO*sw+cO*cw*ci , -cO*si ; ...
      sw*si          ,  cw*si          ,  ci    ];
end

function D = drot_dOm(Om, i, w)
cO=cos(Om); sO=sin(Om); ci=cos(i); si=sin(i); cw=cos(w); sw=sin(w);
D = [ -sO*cw-cO*sw*ci ,  sO*sw-cO*cw*ci , cO*si ; ...
       cO*cw-sO*sw*ci , -cO*sw-sO*cw*ci , sO*si ; ...
       0              ,  0              , 0     ];
end

function D = drot_di(Om, i, w)
cO=cos(Om); sO=sin(Om); ci=cos(i); si=sin(i); cw=cos(w); sw=sin(w);
D = [  sO*sw*si ,  sO*cw*si ,  sO*ci ; ...
      -cO*sw*si , -cO*cw*si , -cO*ci ; ...
       sw*ci    ,  cw*ci    , -si    ];
end

function D = drot_dw(Om, i, w)
cO=cos(Om); sO=sin(Om); ci=cos(i); si=sin(i); cw=cos(w); sw=sin(w);
D = [ -cO*sw-sO*cw*ci , -cO*cw+sO*sw*ci , 0 ; ...
      -sO*sw+cO*cw*ci , -sO*cw-cO*sw*ci , 0 ; ...
       cw*si          , -sw*si          , 0 ];
end
