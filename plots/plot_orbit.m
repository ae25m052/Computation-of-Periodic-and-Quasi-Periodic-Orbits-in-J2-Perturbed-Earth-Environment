function plot_orbit(t,X,T,ttl)
% 3D orbit in canonical units: textured Earth (radius 1 DU), trace
% coloured by time in days.

c = j2_constants();

figure; hold on
he = plot_earth(1);
hg = plot_orbit_gradient(X(:,1:3).', t*c.TU/86400, 'turbo', 1.5);
set([he hg],'HandleVisibility','off');

cb = colorbar;
cb.Label.String = 'time  [days]';

k = t <= T;
plot3(X(k,1),X(k,2),X(k,3),'k','LineWidth',1.5,'DisplayName','first revolution');
legend('Location','northeast');

lim = 1.1*max(sqrt(sum(X(:,1:3).^2,2)));
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);
grid on; axis vis3d; view(45,25)
xlabel('x  [DU]'); ylabel('y  [DU]'); zlabel('z  [DU]');
title(ttl);

end
