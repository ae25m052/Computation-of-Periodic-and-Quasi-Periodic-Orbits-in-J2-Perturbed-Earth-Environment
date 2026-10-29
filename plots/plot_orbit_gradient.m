function h = plot_orbit_gradient(X, t, cmapname, lw)
% Orbit trace coloured by time. X is 3xN (or 6xN), t is 1xN.

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
