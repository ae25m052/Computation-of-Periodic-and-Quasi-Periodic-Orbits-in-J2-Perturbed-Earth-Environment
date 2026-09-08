function [u1, DP, tau, xe] = meridional_return_map(u0, mu, epsilon, J2, Re, hz, E, Tguess, opts)
%MERIDIONAL_RETURN_MAP  First return to the section {z = 0, dz/dt > 0}.
%   [u1, DP, tau, xe] = meridional_return_map(u0, mu, epsilon, J2, Re, hz, E, Tguess, opts)
%       E = 1/2 ( p_rho^2 + p_z^2 + hz^2/rho^2 ) - U(rho,0)
%     =>  p_z = sqrt( 2(E + U) - p_rho^2 - hz^2/rho^2 )     (taking p_z > 0)
%   DP is the projected Jacobian
%       DP = S ( I - f n' / (n' f) ) Phi(tau) Q

rho = u0(1);
pr  = u0(2);

C  = 0.5 * epsilon * J2 * mu * Re^2;
U  = mu/rho + C/rho^3;                       
Up = -mu/rho^2 - 3*C/rho^4;                  

val = 2*(E + U) - pr^2 - hz^2/rho^2;
if val <= 0
    error('meridional_return_map:offSurface', ...
          'No real p_z at rho = %.6f, p_rho = %.6f (energy surface not reached).', rho, pr);
end
pz = sqrt(val);

x0 = [rho ; 0 ; pr ; pz];

Q = zeros(4,2);
Q(1,1) = 1;
Q(3,2) = 1;
Q(4,1) = (Up + hz^2/rho^3) / pz;
Q(4,2) = -pr / pz;

Y0 = [x0 ; reshape(eye(4), 16, 1)];

o = odeset(opts, 'Events', @zsection);
[~, ~, te, Ye, ~] = ode89(@(t,Y) meridional_stm_eom(t, Y, mu, epsilon, J2, Re, hz), ...
                          [0 3*Tguess], Y0, o);

k = find(te > 1e-8, 1, 'first');
if isempty(k)
    error('meridional_return_map:noReturn', 'No upward z = 0 crossing found.');
end
tau = te(k);
Ye  = Ye(k,:).';

xe  = Ye(1:4);
Phi = reshape(Ye(5:20), 4, 4);

u1 = [xe(1) ; xe(3)];

if nargout > 1
    f = meridional_eom(0, xe, mu, epsilon, J2, Re, hz);
    n = [0;1;0;0];
    M = (eye(4) - (f*n.')/(n.'*f)) * Phi;
    S = zeros(2,4); S(1,1) = 1; S(2,3) = 1;
    DP = S * M * Q;
end
end

% -------------------------------------------------------------------------
function [value, isterminal, direction] = zsection(~, Y)
value      = Y(2);      
isterminal = 0;
direction  = 1;         
end
