function dX = rf_eom(t,X,mu,Re,J2,Wf)
% Gravity from j2_eom, plus Coriolis and centrifugal terms.

dX = j2_eom(t,X,mu,Re,J2,1);

dX(4) = dX(4) + 2*Wf*X(5) + Wf^2*X(1);
dX(5) = dX(5) - 2*Wf*X(4) + Wf^2*X(2);

end
