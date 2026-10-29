function h = plot_earth(Re_plot)
% Textured Earth of radius Re_plot (1 in canonical units). The topography
% is converted to RGB so it does not use the axes colormap.

    try
        S = load('topo.mat', 'topo', 'topomap1');
        topo = S.topo;
        cmap = S.topomap1;

        nc   = size(cmap, 1);
        cmin = min(topo(:));
        cmax = max(topo(:));
        idx  = round((topo - cmin) / (cmax - cmin) * (nc - 1)) + 1;
        idx  = min(max(idx, 1), nc);
        rgb  = ind2rgb(idx, cmap);

        [xx, yy, zz] = sphere(80);
        h = surface(Re_plot*xx, Re_plot*yy, Re_plot*zz, ...
                    'FaceColor',    'texturemap', ...
                    'CData',        rgb, ...
                    'EdgeColor',    'none', ...
                    'FaceLighting', 'gouraud');
    catch
        % plain sphere if topo.mat is not available
        [xx, yy, zz] = sphere(60);
        h = surf(Re_plot*xx, Re_plot*yy, Re_plot*zz, ...
                 'FaceColor', [0.25 0.45 0.75], 'EdgeColor', 'none', ...
                 'FaceLighting', 'gouraud');
    end

    material([0.6 0.5 0.1 1 0]);
    axis equal;

    light('Position', [1 1 0.6],     'Style', 'infinite', 'Color', [1 1 0.97]);
    light('Position', [-1 -0.5 -0.3],'Style', 'infinite', 'Color', [0.28 0.28 0.33]);
end
