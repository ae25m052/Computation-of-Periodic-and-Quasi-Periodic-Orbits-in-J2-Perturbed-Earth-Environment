% Inclination sweep from 0 to 180 degrees. Measures the secular drift of
% the node and the apsidal line and compares with first order J2 theory.

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
T  = 2*pi*sqrt(a^3/mu);           % TU

nrev = 100;

inc = [0:1:180, 63.434949, 116.565051];
inc = sort(unique(inc));
N   = numel(inc);

fprintf('a = %.6f DU = %.4f km, e = %.2f, T = %.6f TU = %.4f min\n', ...
        a,a*c.DU,e,T,T*c.TU/60);
fprintf('%d inclinations, %d revolutions each\n\n',N,nrev);

dW = zeros(N,1);  dw = zeros(N,1);
dW_th = zeros(N,1);  dw_th = zeros(N,1);

p = a*(1-e^2);
n = sqrt(mu/a^3);

for k = 1:N
    i = deg2rad(inc(k));
    [r0,v0] = kep2rv(a,e,i,0,0,0,mu);

    [t,X] = ode89(@(t,X) j2_eom(t,X,mu,Re,J2,1), ...
                  linspace(0,nrev*T,nrev*20+1),[r0;v0],c.opts);

    m = numel(t);
    RA = zeros(m,1);  AP = zeros(m,1);
    for q = 1:m
        [~,~,~,RA(q),AP(q)] = rv2kep(X(q,1:3).',X(q,4:6).',mu);
    end
    RA = unwrap(RA);  AP = unwrap(AP);

    % At i = 0 and 180 there is no line of nodes
    if abs(inc(k)) < 1e-8 || abs(inc(k)-180) < 1e-8
        dW(k) = NaN;  dw(k) = NaN;
    else
        dW(k) = rad2deg(RA(end)-RA(1));
        dw(k) = rad2deg(AP(end)-AP(1));
    end

    % first order theory over the same span
    dW_th(k) = rad2deg(-1.5*n*J2*(Re/p)^2*cos(i)        *nrev*T);
    dw_th(k) = rad2deg(0.75*n*J2*(Re/p)^2*(5*cos(i)^2-1)*nrev*T);

    if mod(k,20)==0
        fprintf('  i = %6.2f deg   dOmega = %8.3f   dargp = %8.3f\n', ...
                inc(k),dW(k),dw(k));
    end
end

figure
subplot(2,1,1)
plot(inc,dW,'b',inc,dW_th,'r--'); grid on
ylabel('\Delta\Omega [deg]'); legend('numerical','first order','Location','best')
subplot(2,1,2)
plot(inc,dw,'b',inc,dw_th,'r--'); grid on
xlabel('inclination [deg]'); ylabel('\Delta\omega [deg]')
hold on; plot([63.434949 63.434949],ylim,'k:')

% where the apsidal drift changes sign, looking below 90 deg only
q = find(inc(1:end-1) < 90 & dw(1:end-1).*dw(2:end) < 0, 1);
if ~isempty(q)
    ic = interp1(dw(q:q+1),inc(q:q+1),0);
    fprintf('\napsidal drift vanishes at i = %.4f deg\n',ic);
    fprintf('first order theory gives   i = %.4f deg\n',63.434949);
end

fprintf('max |numerical - theory| in dOmega = %.4f deg\n',max(abs(dW-dW_th)));
fprintf('max |numerical - theory| in dargp  = %.4f deg\n',max(abs(dw-dw_th)));
