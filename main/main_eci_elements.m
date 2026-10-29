% Osculating elements at every instant, integrated in ECI (Wf = 0)

clear; clc; close all

c = j2_constants();

% normalise first
mu = 1;
Re = 1;
J2 = c.J2;

hp = 300;                         % km
e  = 0.65;
rp = (c.Re + hp)/c.DU;            % DU
a  = rp/(1-e);
T  = 2*pi*sqrt(a^3/mu);           % Kepler period, TU

nrev = 2;
inc  = [0 15 30 45 60 75 90];     % deg
t    = linspace(0,nrev*T,1501);

fprintf('ECI, hp = %.1f km, e = %.2f, %d periods, T = %.4f min = %.4f h\n\n',hp,e,nrev,T*c.TU/60,T*c.TU/3600);
fprintf('            peak to peak over the run                change, end - start\n');
fprintf('  i [deg]   a [km]     e          i [deg]          w [deg]     RAAN [deg]\n');

for k = 1:numel(inc)
    i = deg2rad(inc(k));

    [r0,v0] = kep2rv(a,e,i,0,0,0,mu);
    [~,X] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1),t,[r0;v0],c.opts);

    m  = numel(t);
    EL = zeros(m,5);              % a, e, i, w, RAAN
    for q = 1:m
        [EL(q,1),EL(q,2),EL(q,3),EL(q,5),EL(q,4)] = rv2kep(X(q,1:3).',X(q,4:6).',mu);
    end

    % angles to (-180,180] then unwrap, so a start at 0 does not show as 360
    EL(:,4) = unwrap(mod(EL(:,4)+pi,2*pi) - pi);
    EL(:,5) = unwrap(mod(EL(:,5)+pi,2*pi) - pi);

    aK = EL(:,1)*c.DU;
    eK = EL(:,2);
    iD = rad2deg(EL(:,3));
    wD = rad2deg(EL(:,4));
    WD = rad2deg(EL(:,5));

    fprintf('  %5.0f    %8.4f   %.3e   %.3e        %+9.5f   %+9.5f\n', ...
            inc(k),max(aK)-min(aK),max(eK)-min(eK),max(iD)-min(iD), ...
            wD(end)-wD(1),WD(end)-WD(1));

    th  = t*c.TU/3600;            % time in hours
    Th  = T*c.TU/3600;            % one period in hours

    % rates at every instant, in deg per revolution
    %  - from the curves: slope of w(t) and RAAN(t)
    %  - from Brouwer, using the osculating a, e, i at that instant
    wdot = gradient(wD,th)*Th;
    Wdot = gradient(WD,th)*Th;
    nO = sqrt(mu./EL(:,1).^3);
    pO = EL(:,1).*(1 - EL(:,2).^2);
    kO = nO*J2.*(Re./pO).^2*T;    % times T: per revolution
    wB = rad2deg( 0.75*kO.*(5*cos(EL(:,3)).^2 - 1));
    WB = rad2deg(-1.5 *kO.*cos(EL(:,3)));

    lab = {'a [km]','e','i [deg]','\omega [deg]','\Omega [deg]', ...
           'd\omega/dt [deg/rev]','d\Omega/dt [deg/rev]'};
    Y   = [aK eK iD wD WD wdot Wdot];
    figure
    for q = 1:7
        subplot(4,2,q)
        plot(th,Y(:,q),'b'); grid on
        if q == 6, hold on; plot(th,wB,'r--'); end
        if q == 7, hold on; plot(th,WB,'r--'); end
        xlim([0 nrev*Th]); xticks((0:nrev)*Th)
        ylabel(lab{q})
        if q >= 6, xlabel('time [h]'), end
    end
    subplot(4,2,1); title(sprintf('ECI, i = %d deg,  T = %.3f h',inc(k),Th))
    subplot(4,2,6); legend('from curve','Brouwer','Location','best')
    fig_format(gcf)

    plot_orbit(t,X,T,sprintf('ECI, i = %d deg, %d periods',inc(k),nrev));
    fig_format(gcf)
end

fprintf('\nAt i = 0 there is no node: RAAN is set to 0 and w is the perigee\n');
fprintf('direction measured from the ECI x-axis.\n');
