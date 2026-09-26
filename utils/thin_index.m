function idx = thin_index(N, maxpts)
%THIN_INDEX  Index vector subsampling 1:N down to at most maxpts entries.
%
%   idx = thin_index(numel(t), 12000);
%   plot3(X(1,idx), X(2,idx), X(3,idx));
%
%   Four hundred revolutions leave tens of thousands of integrator steps.
%   Handing all of them to plot3 -- and worse, to the patch object inside
%   plot_orbit_gradient -- makes MATLAB crawl and bloats the exported PNG
%   without adding anything the eye can resolve at 300 DPI.
%
%   The first and last points are always kept, so the trace still starts
%   and ends exactly where the trajectory does.
%
%   M.Tech. project, Raj Khismatrao (AE25M052).

    if N <= maxpts
        idx = 1:N;
    else
        idx = unique([1, round(linspace(1, N, maxpts)), N]);
    end
end
