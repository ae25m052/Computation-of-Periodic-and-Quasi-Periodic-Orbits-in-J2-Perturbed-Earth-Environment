function fig_format(fig,w)
% Report formatting for a figure: white background, bold 14 pt text,
% line width 1.5. w = figure width in cm (default 16, about the report
% text width, so the fonts are not shrunk when the figure is inserted).

if nargin < 1, fig = gcf; end
if nargin < 2, w = 16;   end

theme(fig,'light');                  % white figure even in dark mode

fig.WindowStyle = 'normal';          % size cannot be set while docked
fig.Units       = 'centimeters';
fig.Position    = [2 2 w 0.75*w];

ax = findobj(fig,'Type','axes');
set(ax,'FontSize',14,'FontWeight','bold','LineWidth',1.5,'GridAlpha',0.25)
set([ax.XLabel ax.YLabel ax.ZLabel ax.Title],'FontWeight','bold')

set(findobj(fig,'Type','line'),'LineWidth',1.5)
set(findobj(fig,'Type','constantline'),'LineWidth',1.5)
set(findobj(fig,'Type','patch'),'LineWidth',1.5)

set(findobj(fig,'Type','legend'),'FontSize',12,'FontWeight','bold')
set(findobj(fig,'Type','colorbar'),'FontSize',12,'FontWeight','bold','LineWidth',1.5)

end
