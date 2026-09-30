function A = j2_A_matrix(r,mu,Re,J2,eps)

x = r(1);  y = r(2);  z = r(3);
rn = norm(r);

% two-body part
d = -mu/rn^3;
c =  3*mu/rn^5;

Gxx = d + c*x*x;  Gyy = d + c*y*y;  Gzz = d + c*z*z;
Gxy =     c*x*y;  Gxz =     c*x*z;  Gyz =     c*y*z;

% J2 part
k  = 1.5*J2*mu*Re^2;
f  = 1/rn^5 - 5*z^2/rn^7;
g  = 3/rn^5 - 5*z^2/rn^7;
fx = -5*x/rn^7  + 35*z^2*x/rn^9;
fy = -5*y/rn^7  + 35*z^2*y/rn^9;
fz = -15*z/rn^7 + 35*z^3/rn^9;
gz = -25*z/rn^7 + 35*z^3/rn^9;

Gxx = Gxx - eps*k*(f + x*fx);
Gyy = Gyy - eps*k*(f + y*fy);
Gzz = Gzz - eps*k*(g + z*gz);
Gxy = Gxy - eps*k*x*fy;
Gxz = Gxz - eps*k*x*fz;
Gyz = Gyz - eps*k*y*fz;

G = [Gxx Gxy Gxz
     Gxy Gyy Gyz
     Gxz Gyz Gzz];

A = [zeros(3) eye(3); G zeros(3)];

end
