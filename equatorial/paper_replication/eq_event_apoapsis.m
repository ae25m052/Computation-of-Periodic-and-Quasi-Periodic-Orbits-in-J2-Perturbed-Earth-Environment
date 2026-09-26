function [value, isterminal, direction] = eq_event_apoapsis(~, Y)
%EQ_EVENT_APOAPSIS  ode89 event: the orbit reaches its farthest point.
%
%   The radial velocity is  rdot = (r . v) / |r|.  Its sign equals the sign
%   of r . v = x*xdot + y*ydot + z*zdot, so we watch that instead.
%
%   Starting at periapsis, r . v = 0 and then grows positive (moving
%   outward). The first time it falls back through zero, going from + to
%   -, the orbit is at apoapsis. That is HALF a radial period.
%
%   Used by eq_dc_pseudo.m and main_equatorial_dc.m.

value      = Y(1)*Y(4) + Y(2)*Y(5) + Y(3)*Y(6);   % r . v
isterminal = 1;
direction  = -1;
end
