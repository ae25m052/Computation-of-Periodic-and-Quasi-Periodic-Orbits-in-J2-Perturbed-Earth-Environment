function dX = j2_eom(~,X,mu,Re,J2,eps)
% Two-body plus eps*J2 acceleration. eps = 0 gives Kepler, eps = 1 full J2.

r = X(1:3);
v = X(4:6);

rn = norm(r);
z  = r(3);

k   = 1.5*J2*mu*Re^2/rn^5;
aJ2 = -k*[ r(1)*(1 - 5*z^2/rn^2)
           r(2)*(1 - 5*z^2/rn^2)
           z   *(3 - 5*z^2/rn^2) ];

dX = [v; -mu*r/rn^3 + eps*aJ2];

end
