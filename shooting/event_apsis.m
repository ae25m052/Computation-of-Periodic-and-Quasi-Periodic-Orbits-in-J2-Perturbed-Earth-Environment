function [value,isterminal,direction] = event_apsis(~,X)
% Event at an apsis. r.v is zero there

value      = X(1:3).'*X(4:6);
isterminal = 1;
direction  = -1;

end
