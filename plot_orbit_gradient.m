function h = plot_orbit_gradient(X, t, cmapname, lw)
%   Inputs:
%       X        - 3xN or 6xN state history (rows 1:3 used)
%       t        - 1xN time vector, same length as columns of X, used
%                   only to color the trace early->late
%       cmapname - colormap name, e.g. 'turbo', 'parula', 'cool'
%       lw       - line width
%
%   Output:
%       h - handle to the patch object

    if nargin < 3, cmapname = 'turbo'; end
    if nargin < 4, lw = 1.5; end

    x = X(1,:); y = X(2,:); z = X(3,:);

    h = patch([x NaN], [y NaN], [z NaN], [t(:).' NaN], ...
              'FaceColor', 'none', 'EdgeColor', 'interp', ...
              'LineWidth', lw, 'Marker', 'none');
    colormap(gca, cmapname);

    if numel(t) > 1 && (max(t) > min(t))
        if exist('clim', 'file') || exist('clim', 'builtin')
            clim(gca, [min(t) max(t)]);
        else
            caxis(gca, [min(t) max(t)]); 
        end
    end
end
