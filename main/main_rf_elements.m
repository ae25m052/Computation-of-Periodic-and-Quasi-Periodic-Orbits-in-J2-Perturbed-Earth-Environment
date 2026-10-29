% Osculating elements at every instant in the rotating frame

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
p  = a*(1-e^2);
n  = sqrt(mu/a^3);
T  = 2*pi*sqrt(a^3/mu);           % Kepler period, TU

nrev = 2;
inc  = [0 15 30 45 60 75 90];     % deg
t    = linspace(0,nrev*T,1501);

fprintf('hp = %.1f km, e = %.2f, %d periods, T = %.4f min = %.4f h\n\n',hp,e,nrev,T*c.TU/60,T*c.TU/3600);
fprintf('                         peak to peak over the run                change, end - start\n');
fprintf('  i [deg]  Wf [deg/day]   a [km]     e          i [deg]          w [deg]     RAAN [deg]\n');

for k = 1:numel(inc)
    i  = deg2rad(inc(k));
    Wf = -1.5*n*J2*(Re/p)^2*cos(i);
    W  = [0; 0; Wf];

    [r0,v0] = kep2rv(a,e,i,0,0,0,mu);
    X0 = [r0; v0 - cross(W,r0)];
    [~,X] = ode89(@(t,X) rf_eom(t,X,mu,Re,J2,Wf),t,X0,c.opts);

    m  = numel(t);
    EL = zeros(m,5);              % a, e, i, w, RAAN
    for q = 1:m
        r = X(q,1:3).';
        v = X(q,4:6).' + cross(W,r);      % inertial velocity
        [EL(q,1),EL(q,2),EL(q,3),EL(q,5),EL(q,4)] = rv2kep(r,v,mu);
    end

    EL(:,4) = unwrap(mod(EL(:,4)+pi,2*pi) - pi);
    EL(:,5) = unwrap(mod(EL(:,5)+pi,2*pi) - pi);

    aK = EL(:,1)*c.DU;
    eK = EL(:,2);
    iD = rad2deg(EL(:,3));
    wD = rad2deg(EL(:,4));
    WD = rad2deg(EL(:,5));

    fprintf('  %5.0f    %+9.4f    %8.4f   %.3e   %.3e        %+9.5f   %+9.5f\n', ...
            inc(k),rad2deg(Wf/c.TU*86400),max(aK)-min(aK),max(eK)-min(eK), ...
            max(iD)-min(iD),wD(end)-wD(1),WD(end)-WD(1));

    th = t*c.TU/3600;             % time in hours
    Th = T*c.TU/3600;             % one period in hours
    lab = {'a [km]','e','i [deg]','\omega [deg]','\Omega [deg]'};
    Y   = [aK eK iD wD WD];
    figure
    for q = 1:5
        subplot(3,2,q)
        plot(th,Y(:,q),'b'); grid on
        xlim([0 nrev*Th]); xticks((0:nrev)*Th)
        ylabel(lab{q})
        if q >= 4, xlabel('time [h]'), end
    end
    subplot(3,2,1); title(sprintf('i = %d deg,  T = %.3f h',inc(k),Th))
    fig_format(gcf)

    % the orbit itself, rotating frame, same 2 periods
    plot_orbit(t,X,T,sprintf('Rotating frame, i = %d deg, %d periods',inc(k),nrev));
    fig_format(gcf)

    % radius difference from the Kepler orbit with the same start
    [~,XK] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,0),t,[r0;v0],c.opts);
    dr = (sqrt(sum(X(:,1:3).^2,2)) - sqrt(sum(XK(:,1:3).^2,2)))*c.DU;

    figure
    plot(th,dr,'b'); grid on
    xlim([0 nrev*Th]); xticks((0:nrev)*Th)
    xlabel('time [h]'); ylabel('r_{J2} - r_{Kepler} [km]')
    title(sprintf('i = %d deg',inc(k)))
    fig_format(gcf)
end

fprintf('\nAt i = 0 there is no node: RAAN is set to 0 and w is the perigee\n');
fprintf('direction measured from the rotating x-axis.\n');
