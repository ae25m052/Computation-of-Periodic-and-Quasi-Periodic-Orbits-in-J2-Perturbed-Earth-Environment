function [a,e,i,RAAN,argp,nu] = rv2kep(r,v,mu)
% Position and velocity to classical elements. Angles in radians.

rmag = norm(r);
vmag = norm(v);

h  = cross(r,v);   hmag = norm(h);
n  = cross([0;0;1],h);   nmag = norm(n);

evec = ((vmag^2 - mu/rmag)*r - dot(r,v)*v)/mu;
e    = norm(evec);

a = 1/(2/rmag - vmag^2/mu);
i = acos(h(3)/hmag);

if nmag > 1e-10
    RAAN = atan2(n(2),n(1));
    argp = atan2(dot(cross(n,evec),h)/hmag, dot(n,evec));
else
    RAAN = 0;
    argp = atan2(evec(2),evec(1));
end

nu = atan2(dot(cross(evec,r),h)/hmag, dot(evec,r));

RAAN = mod(RAAN,2*pi);
argp = mod(argp,2*pi);
nu   = mod(nu,2*pi);

end
