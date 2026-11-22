function [u,V,Wf,T,it,DP,cT] = shoot_circular(r0,i,u,V,Wf,mu,Re,J2,opts)
% Near-circular periodic orbit
% DP = Poincare-map matrix (from the monodromy matrix), cT = lap-time row
it = [];
for k = 0:20
    X0 = [r0; 0; 0; u; V*cos(i) - Wf*r0; V*sin(i)];
    [XS,T,DP,psiS,cT] = one_lap(X0,mu,Re,J2,Wf,opts);

    % miss: x, y and xdot must come back (ydot, zdot follow from C and Hz)
    F  = XS([1 2 4]) - X0([1 2 4]);
    it = [it; k, u, V, Wf, T, XS.', norm(F)];
    if norm(F) < 1e-12, return; end

    % how the start state changes with each unknown
    dXu = [0; 0; 0; 1; 0; 0];
    dXV = [0; 0; 0; 0; cos(i); sin(i)];
    dXW = [0; 0; 0; 0; -r0; 0];

    % columns of J: crossing change minus start change
    J  = [DP*dXu - dXu, DP*dXV - dXV, DP*dXW + psiS - dXW];
    dq = J([1 2 4],:)\F;

    u  = u  - dq(1);
    V  = V  - dq(2);
    Wf = Wf - dq(3);
end
warning('shoot_circular: no convergence in 20 steps')

end


function [XS,T,DP,psiS,cT] = one_lap(X0,mu,Re,J2,Wf,opts)
% From the node to the next upward crossing of z = 0, with Phi and psi.

Y0 = [X0; reshape(eye(6),36,1); zeros(6,1)];
Tc = 2*pi*sqrt(norm(X0(1:3))^3/mu);      % circular period, only for time spans
f  = @(t,Y) rf_stm_eom(t,Y,mu,Re,J2,Wf);

% first piece without events, so the start crossing is not counted
[t1,Y1] = ode89(f,[0 0.25*Tc],Y0,opts);
o = odeset(opts,'Events',@event_node);
[~,~,te,Ye] = ode89(f,[t1(end) 10*Tc],Y1(end,:).',o);

T   = te(1);
Y   = Ye(1,:).';
XS  = Y(1:6);
Phi = reshape(Y(7:42),6,6);              % monodromy matrix at the fixed point
psi = Y(43:48);

% remove the shift in crossing time
Xd   = rf_eom(T,XS,mu,Re,J2,Wf);
DP   = Phi - Xd*Phi(3,:)/Xd(3);
psiS = psi - Xd*psi(3)/Xd(3);
cT   = -Phi(3,:)/Xd(3);                  % lap-time change: dT = cT*dX

end


function [val,stop,dir] = event_node(~,Y)
val  = Y(3);     % z = 0
stop = 1;
dir  = 1;        % going up only
end
