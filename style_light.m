function style_light(figHandle, fontSize)
%   Inputs:
%       figHandle - figure handle (default: gcf)
%       fontSize  - base font size in points (default: 12)

    if nargin < 1 || isempty(figHandle), figHandle = gcf;    end
    if nargin < 2 || isempty(fontSize),  fontSize  = 12;     end

    dark = [0.10 0.10 0.10];

    set(figHandle, 'Color', 'w', 'InvertHardcopy', 'off');

    % ---- axes ----
    ax = findall(figHandle, 'Type', 'axes');
    for k = 1:numel(ax)
        a = ax(k);
        set(a, 'Color',     'w', ...
               'XColor',    dark, 'YColor', dark, 'ZColor', dark, ...
               'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.18, ...
               'MinorGridColor', [0.3 0.3 0.3], ...
               'FontName',  'Calibri', ...
               'FontSize',  fontSize, ...
               'LineWidth', 1.0);
        set([a.XLabel a.YLabel a.ZLabel], ...
            'Color', dark, 'FontName', 'Calibri', 'FontSize', fontSize);
        set(a.Title, 'Color', dark, 'FontName', 'Calibri', ...
                     'FontSize', fontSize + 1, 'FontWeight', 'bold');
    end

    % ---- legends ----
    lg = findall(figHandle, 'Type', 'legend');
    for k = 1:numel(lg)
        set(lg(k), 'TextColor', dark, 'Color', 'w', ...
                   'EdgeColor', [0.6 0.6 0.6], ...
                   'FontName', 'Calibri', 'FontSize', max(fontSize-1, 9));
    end

    % ---- colorbars ----
    cb = findall(figHandle, 'Type', 'colorbar');
    for k = 1:numel(cb)
        set(cb(k), 'Color', dark, 'FontName', 'Calibri', ...
                   'FontSize', max(fontSize-1, 9));
        if ~isempty(cb(k).Label)
            set(cb(k).Label, 'Color', dark, 'FontName', 'Calibri');
        end
    end

    allObj = findall(figHandle);
    for k = 1:numel(allObj)
        o = allObj(k);
        if isa(o, 'matlab.graphics.axis.Axes') || ...
           isa(o, 'matlab.graphics.illustration.Legend') || ...
           isa(o, 'matlab.graphics.illustration.ColorBar')
            continue
        end
        if isprop(o, 'String') && isprop(o, 'Color')
            cls = class(o);
            if contains(cls, 'subplot.Text')

                set(o, 'Color', dark, 'FontName', 'Calibri', ...
                       'FontSize', fontSize + 3, 'FontWeight', 'bold');
            else
                c = get(o, 'Color');
                if isnumeric(c) && numel(c) == 3 && all(c > 0.45)
                    set(o, 'Color', dark);   % was light-on-dark; flip it
                end
            end
        end
    end

    drawnow;
end
