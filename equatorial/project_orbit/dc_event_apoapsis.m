function [value, isterminal, direction] = dc_event_apoapsis(~, Y)
%DC_EVENT_APOAPSIS  ode89 event: the orbit reaches apoapsis (farthest point).
%
%   Radial velocity rdot has the sign of r . v = x*xdot + y*ydot + z*zdot.
%   Starting at perigee r . v = 0 and then grows positive (moving outward).
%   The first time it falls back through zero, from + to -, the satellite
%   is at apoapsis: half a radial period after perigee.

value      = Y(1)*Y(4) + Y(2)*Y(5) + Y(3)*Y(6);
isterminal = 1;          % stop there
direction  = -1;         % only + to -
end
