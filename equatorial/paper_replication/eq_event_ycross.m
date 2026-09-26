function [value, isterminal, direction] = eq_event_ycross(~, Y)
%EQ_EVENT_YCROSS  ode89 event: the orbit crosses the x-axis going downward.
%
%   The orbit starts on the +x axis (y = 0) moving in +y. The first time y
%   comes back to zero while DEcreasing is the crossing of the -x axis,
%   i.e. HALF a revolution in angle (theta = 180 deg).
%
%   Used by eq_dc_circular.m.
%
%   Y is the 42-element state+STM vector; only Y(2) = y is used.

value      = Y(2);     % y
isterminal = 1;        % stop the integration here
direction  = -1;       % only when y goes from + to -
end
