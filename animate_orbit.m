function animate_orbit(X, Re_plot, titleStr, gifname, skip)

    pos = X(1:3, 1:skip:end);
    N = size(pos, 2);

    fig = figure('Name', titleStr, 'Color', 'w');
    plot_earth(Re_plot);
    hold on;
    h_line = animatedline('Color', [0.85 0.15 0.15], 'LineWidth', 1.6);
    h_pt   = plot3(pos(1,1), pos(2,1), pos(3,1), 'o', ...
                    'MarkerFaceColor', [1 0.8 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
    grid on; set(gca, 'GridAlpha', 0.25);
    lim = max(abs(pos(:))) * 1.1;
    xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);
    xlabel('x'); ylabel('y'); zlabel('z');
    title(titleStr, 'FontWeight', 'bold');
    view(45, 25); axis vis3d;

    darkTxt = [0.10 0.10 0.10];
    set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
    ax = gca;
    set(ax, 'Color', 'w', ...
            'XColor', darkTxt, 'YColor', darkTxt, 'ZColor', darkTxt, ...
            'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.18, ...
            'FontName', 'Calibri', 'FontSize', 12);
    set([ax.XLabel ax.YLabel ax.ZLabel], 'Color', darkTxt, 'FontName', 'Calibri');
    set(ax.Title, 'Color', darkTxt, 'FontName', 'Calibri', ...
                  'FontSize', 14, 'FontWeight', 'bold');

    if exist('style_light', 'file')
        style_light(fig, 12);
        set(ax.Title, 'Color', darkTxt);   % re-assert after styling
    end

    for k = 1:N
        addpoints(h_line, pos(1,k), pos(2,k), pos(3,k));
        set(h_pt, 'XData', pos(1,k), 'YData', pos(2,k), 'ZData', pos(3,k));
        drawnow limitrate;

        frame = getframe(fig);
        img = frame2im(frame);
        [imind, cm] = rgb2ind(img, 256);
        if k == 1
            imwrite(imind, cm, gifname, 'gif', 'Loopcount', inf, 'DelayTime', 0.03);
        else
            imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime', 0.03);
        end
    end

    fprintf('Saved animation: %s\n', gifname);
end
