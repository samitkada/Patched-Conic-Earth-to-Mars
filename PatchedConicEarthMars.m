%% ============================================================
%  PATCHED CONIC INTERPLANETARY TRAJECTORY: EARTH TO MARS
%  Three-phase model:
%    Phase 1 — Hyperbolic escape from Earth
%    Phase 2 — Heliocentric transfer ellipse (Lambert / Hohmann)
%    Phase 3 — Hyperbolic capture at Mars
%
%  The "patched conic" method stitches these three conics together
%  at the sphere of influence (SOI) boundaries of each planet.
%  Inside a planet's SOI, only that planet's gravity matters.
%  Outside both SOIs, only the Sun's gravity matters.
%
%  Tools: MATLAB | Author: [Your Name]
%% ============================================================

clc; clear; close all;

%% ============================================================
%  SECTION 1: CONSTANTS & PLANETARY DATA
%% ============================================================

% Gravitational parameters (km^3/s^2)
mu_sun   = 1.327e11;    % Sun
mu_earth = 3.986e5;     % Earth
mu_mars  = 4.283e4;     % Mars

% Planet radii (km)
R_earth  = 6378;
R_mars   = 3396;

% Semi-major axes of planetary orbits (km) — circular approximation
a_earth  = 1.496e8;     % 1 AU
a_mars   = 2.279e8;     % 1.524 AU

% Sphere of influence radii (km)
% SOI = a_planet * (m_planet/m_sun)^(2/5)
SOI_earth = 9.25e5;     % ~925,000 km
SOI_mars  = 5.77e5;     % ~577,000 km

% Circular orbital velocities (km/s)
v_earth  = sqrt(mu_sun / a_earth);   % ~29.78 km/s
v_mars   = sqrt(mu_sun / a_mars);    % ~24.13 km/s

% Parking orbit altitudes (km above surface)
h_park_earth = 200;
h_park_mars  = 300;

r_park_earth = R_earth + h_park_earth;
r_park_mars  = R_mars  + h_park_mars;

% Parking orbit velocities (km/s)
v_park_earth = sqrt(mu_earth / r_park_earth);
v_park_mars  = sqrt(mu_mars  / r_park_mars);

%% ============================================================
%  SECTION 2: HELIOCENTRIC TRANSFER ELLIPSE (Hohmann)
%  Minimum-energy transfer between Earth and Mars orbits
%  This is the spine of the patched conic trajectory
%% ============================================================

a_transfer   = (a_earth + a_mars) / 2;  % semi-major axis of transfer ellipse
e_transfer   = (a_mars - a_earth) / (a_mars + a_earth);  % eccentricity

% Velocities at perihelion (Earth) and aphelion (Mars) of transfer ellipse
v_trans_peri = sqrt(mu_sun * (2/a_earth - 1/a_transfer));  % at Earth SOI exit
v_trans_aph  = sqrt(mu_sun * (2/a_mars  - 1/a_transfer));  % at Mars SOI entry

% Heliocentric velocity excess at each planet
v_inf_earth  = v_trans_peri - v_earth;  % hyperbolic excess at Earth (km/s)
v_inf_mars   = v_mars - v_trans_aph;    % hyperbolic excess at Mars  (km/s)

% Transfer time (half-period of transfer ellipse)
T_transfer   = pi * sqrt(a_transfer^3 / mu_sun);  % seconds
T_days       = T_transfer / 86400;                % days

%% ============================================================
%  SECTION 3: PHASE 1 — HYPERBOLIC ESCAPE FROM EARTH
%  The spacecraft leaves its parking orbit with a delta-V burn,
%  achieving v_inf relative to Earth, and exits Earth's SOI
%  on a hyperbolic trajectory
%% ============================================================

% Hyperbolic excess speed relative to Earth
C3_earth     = v_inf_earth^2;   % characteristic energy (km^2/s^2)

% Velocity at periapsis of departure hyperbola (at parking orbit radius)
v_hyp_earth  = sqrt(C3_earth + 2*mu_earth/r_park_earth);

% Delta-V for Trans-Mars Injection (TMI) burn
dv_TMI       = v_hyp_earth - v_park_earth;

% Hyperbola parameters at Earth
a_hyp_earth  = -mu_earth / C3_earth;           % negative (hyperbola)
e_hyp_earth  =  1 + r_park_earth / abs(a_hyp_earth);

% Generate departure hyperbola trajectory
% True anomaly range: from periapsis to SOI boundary
nu_SOI_earth = acos((a_hyp_earth*(1-e_hyp_earth^2)/SOI_earth - 1) / e_hyp_earth);
nu_e         = linspace(0, nu_SOI_earth * 0.92, 300);
r_hyp_e      = a_hyp_earth*(1-e_hyp_earth^2) ./ (1 + e_hyp_earth*cos(nu_e));
r_hyp_e      = abs(r_hyp_e);

x_hyp_e      = r_hyp_e .* cos(nu_e);
y_hyp_e      = r_hyp_e .* sin(nu_e);

%% ============================================================
%  SECTION 4: PHASE 3 — HYPERBOLIC CAPTURE AT MARS
%  Spacecraft arrives at Mars SOI with v_inf, fires retro burn
%  to slow into a capture orbit (Mars Orbit Insertion, MOI)
%% ============================================================

C3_mars      = v_inf_mars^2;
v_hyp_mars   = sqrt(C3_mars + 2*mu_mars/r_park_mars);
dv_MOI       = v_hyp_mars - v_park_mars;

% Hyperbola parameters at Mars
a_hyp_mars   = -mu_mars / C3_mars;
e_hyp_mars   =  1 + r_park_mars / abs(a_hyp_mars);

% Generate arrival hyperbola (approach from SOI to periapsis)
nu_SOI_mars  = acos((a_hyp_mars*(1-e_hyp_mars^2)/SOI_mars - 1) / e_hyp_mars);
nu_m         = linspace(nu_SOI_mars * 0.92, 0, 300);
r_hyp_m      = a_hyp_mars*(1-e_hyp_mars^2) ./ (1 + e_hyp_mars*cos(nu_m));
r_hyp_m      = abs(r_hyp_m);

x_hyp_m      = r_hyp_m .* cos(nu_m);
y_hyp_m      = r_hyp_m .* sin(nu_m);

%% ============================================================
%  SECTION 5: HELIOCENTRIC TRAJECTORIES
%  Numerical propagation of transfer ellipse using ode45
%  Two-body equations: d²r/dt² = -mu/|r|³ * r
%% ============================================================

% Initial state on transfer ellipse at Earth's position
% Spacecraft at perihelion: position = [a_earth, 0], velocity = [0, v_trans_peri]
state0 = [a_earth; 0; 0; v_trans_peri];   % [x, y, vx, vy] in km, km/s

ode_2body = @(t, s) [s(3); s(4);
                     -mu_sun*s(1)/(s(1)^2+s(2)^2)^1.5;
                     -mu_sun*s(2)/(s(1)^2+s(2)^2)^1.5];

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t_traj, S_traj] = ode45(ode_2body, [0 T_transfer], state0, opts);

x_trans = S_traj(:,1);
y_trans = S_traj(:,2);

% Earth's full orbit
theta_full = linspace(0, 2*pi, 500);
x_earth_orb = a_earth * cos(theta_full);
y_earth_orb = a_earth * sin(theta_full);

% Mars's full orbit
x_mars_orb  = a_mars * cos(theta_full);
y_mars_orb  = a_mars * sin(theta_full);

% Mars position at arrival (rotated by pi from Earth departure)
mars_arrival_angle = pi;   % at aphelion of transfer ellipse
x_mars_arr = a_mars * cos(mars_arrival_angle);
y_mars_arr = a_mars * sin(mars_arrival_angle);

%% ============================================================
%  SECTION 6: RESULTS
%% ============================================================

total_dv = abs(dv_TMI) + abs(dv_MOI);

fprintf('\n=====================================================\n');
fprintf('   PATCHED CONIC: EARTH TO MARS TRAJECTORY\n');
fprintf('=====================================================\n');
fprintf('HELIOCENTRIC TRANSFER\n');
fprintf('  Transfer semi-major axis:  %.3e km\n',  a_transfer);
fprintf('  Transfer eccentricity:     %.4f\n',      e_transfer);
fprintf('  Transfer time:             %.1f days\n', T_days);
fprintf('  v at Earth (perihelion):   %.4f km/s\n', v_trans_peri);
fprintf('  v at Mars  (aphelion):     %.4f km/s\n', v_trans_aph);
fprintf('-----------------------------------------------------\n');
fprintf('EARTH DEPARTURE (Phase 1)\n');
fprintf('  Parking orbit altitude:    %d km\n',     h_park_earth);
fprintf('  Parking orbit velocity:    %.4f km/s\n', v_park_earth);
fprintf('  v_inf (hyperbolic excess): %.4f km/s\n', v_inf_earth);
fprintf('  C3 (characteristic nrgy): %.4f km2/s2\n',C3_earth);
fprintf('  Delta-V (TMI burn):        %.4f km/s\n', dv_TMI);
fprintf('-----------------------------------------------------\n');
fprintf('MARS ARRIVAL (Phase 3)\n');
fprintf('  Capture orbit altitude:    %d km\n',     h_park_mars);
fprintf('  v_inf at Mars SOI:         %.4f km/s\n', v_inf_mars);
fprintf('  Delta-V (MOI burn):        %.4f km/s\n', dv_MOI);
fprintf('-----------------------------------------------------\n');
fprintf('TOTAL MISSION Delta-V:       %.4f km/s\n', total_dv);
fprintf('=====================================================\n\n');

%% ============================================================
%  SECTION 7: 3D PLOTS — SINGLE FIGURE
%% ============================================================

ax_bg  = [0.07 0.07 0.12];
col_e  = [0.3  0.6  1.0];
col_m  = [1.0  0.45 0.2];
col_sun= [1.0  0.85 0.2];
col_hy = [0.4  1.0  0.6];
col_bv = [1.0  1.0  0.2];

theta_circ = linspace(0, 2*pi, 200);
inc_mars   = deg2rad(1.85);

sphere_fn = @(R, xc, yc, zc) deal(...
    xc + R*sin(linspace(0,pi,40)') * cos(linspace(0,2*pi,40)), ...
    yc + R*sin(linspace(0,pi,40)') * sin(linspace(0,2*pi,40)), ...
    zc + R*cos(linspace(0,pi,40)') * ones(1,40));

fig = figure('Name','Patched Conic — Earth to Mars', ...
             'Color',[0.04 0.04 0.08], ...
             'Position',[60 40 1400 820]);

%% ============================================================
%  LEFT (large): Heliocentric 3D view
%% ============================================================
ax1 = subplot(1,3,[1 2]);
set(ax1,'Color',ax_bg,...
    'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5],'ZColor',[0.5 0.5 0.5],...
    'GridColor',[0.15 0.15 0.22],'GridAlpha',1,...
    'XGrid','on','YGrid','on','ZGrid','on',...
    'FontSize',9,'FontName','Consolas');
hold on; axis equal;

% Sun
[xs,ys,zs] = sphere_fn(0.04e8, 0, 0, 0);
surf(xs,ys,zs,'FaceColor',col_sun,'EdgeColor','none','FaceAlpha',0.95);
text(0,0,0.07e8,'Sun','Color',col_sun,'FontSize',9,...
     'HorizontalAlignment','center','FontName','Consolas');

% Orbit rings
plot3(x_earth_orb, y_earth_orb, zeros(size(x_earth_orb)),...
      '--','Color',[col_e 0.35],'LineWidth',0.8);
z_mars_orb = a_mars*sin(inc_mars)*sin(theta_full);
plot3(x_mars_orb, y_mars_orb, z_mars_orb,...
      '--','Color',[col_m 0.35],'LineWidth',0.8);

% Transfer trajectory colored by time
n_seg = length(x_trans)-1;
cmap  = cool(n_seg);
z_trans_line = linspace(0, a_mars*sin(inc_mars)*sin(pi), length(x_trans));
for k = 1:n_seg
    plot3(x_trans(k:k+1), y_trans(k:k+1), z_trans_line(k:k+1),...
          '-','Color',cmap(k,:),'LineWidth',2.5);
end

% Earth sphere
[xe,ye,ze] = sphere_fn(0.025e8, a_earth, 0, 0);
surf(xe,ye,ze,'FaceColor',col_e,'EdgeColor','none','FaceAlpha',0.85);
text(a_earth, 0.12e8, 0.05e8,'Earth','Color',col_e,'FontSize',9,'FontName','Consolas');

% Mars sphere
z_mars_a = a_mars*sin(inc_mars)*sin(mars_arrival_angle);
[xm,ym,zm] = sphere_fn(0.018e8, x_mars_arr, y_mars_arr, z_mars_a);
surf(xm,ym,zm,'FaceColor',col_m,'EdgeColor','none','FaceAlpha',0.85);
text(x_mars_arr-0.08e8, y_mars_arr, z_mars_a+0.07e8,'Mars','Color',col_m,...
     'FontSize',9,'FontName','Consolas');

% Burn markers
scatter3(a_earth, 0, 0, 120, col_bv,'filled','MarkerEdgeColor','w','LineWidth',1);
scatter3(x_mars_arr, y_mars_arr, z_mars_a, 120, col_bv,'filled','MarkerEdgeColor','w','LineWidth',1);
text(a_earth+0.05e8,-0.12e8,0.02e8,sprintf('\\DeltaV_1=%.2f km/s',dv_TMI),...
     'Color',col_bv,'FontSize',8,'FontName','Consolas');
text(x_mars_arr-0.22e8,y_mars_arr,z_mars_a+0.05e8,sprintf('\\DeltaV_2=%.2f km/s',dv_MOI),...
     'Color',col_bv,'FontSize',8,'FontName','Consolas');

% Drop lines to ecliptic for depth
plot3([a_earth a_earth],[0 0],[0 0],'--','Color',[col_e 0.25],'LineWidth',0.6);
plot3([x_mars_arr x_mars_arr],[y_mars_arr y_mars_arr],[0 z_mars_a],...
      '--','Color',[col_m 0.25],'LineWidth',0.6);

colormap(ax1, cool);
cb1 = colorbar(ax1);
cb1.Ticks=[0 1]; cb1.TickLabels={'Departure','Arrival'};
cb1.Color=[0.6 0.6 0.6];
cb1.Label.String=sprintf('TOF  (%.0f days)',T_days);
cb1.Label.Color=[0.6 0.6 0.6]; cb1.Label.FontName='Consolas';
xlabel('X (km)','FontName','Consolas');
ylabel('Y (km)','FontName','Consolas');
zlabel('Z (km)','FontName','Consolas');
view(-38, 22);
title('Heliocentric transfer','Color','w','FontSize',11,'FontName','Consolas');

%% ============================================================
%  TOP RIGHT: Earth departure hyperbola
%% ============================================================
ax2 = subplot(2,3,3);
set(ax2,'Color',ax_bg,...
    'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5],'ZColor',[0.5 0.5 0.5],...
    'GridColor',[0.15 0.15 0.22],'GridAlpha',1,...
    'XGrid','on','YGrid','on','ZGrid','on',...
    'FontSize',8,'FontName','Consolas');
hold on; axis equal;

[xe2,ye2,ze2] = sphere_fn(R_earth, 0, 0, 0);
surf(xe2,ye2,ze2,'FaceColor',col_e,'EdgeColor','none','FaceAlpha',0.8);
[xa,ya,za_atm] = sphere_fn(R_earth*1.05, 0, 0, 0);
surf(xa,ya,za_atm,'FaceColor',col_e,'EdgeColor','none','FaceAlpha',0.07);

plot3(r_park_earth*cos(theta_circ), r_park_earth*sin(theta_circ),...
      zeros(size(theta_circ)),'--','Color',[col_e 0.6],'LineWidth',1);
plot3(SOI_earth*cos(theta_circ), SOI_earth*sin(theta_circ),...
      zeros(size(theta_circ)),':','Color',[0.5 0.5 0.5 0.6],'LineWidth',0.7);
text(SOI_earth*0.65,SOI_earth*0.65,SOI_earth*0.05,'SOI',...
     'Color',[0.5 0.5 0.5],'FontSize',7,'FontName','Consolas');

z_hyp_e = linspace(0, SOI_earth*0.18, length(x_hyp_e));
plot3(x_hyp_e, y_hyp_e, z_hyp_e,'Color',col_hy,'LineWidth',2);

scatter3(r_park_earth,0,0,80,col_bv,'filled','MarkerEdgeColor','w');
text(r_park_earth*1.1,r_park_earth*0.5,r_park_earth*0.2,...
     sprintf('TMI  \\DeltaV=%.2f',dv_TMI),...
     'Color',col_bv,'FontSize',7,'FontName','Consolas');

quiver3(x_hyp_e(end),y_hyp_e(end),z_hyp_e(end),...
        cos(nu_SOI_earth*0.85)*SOI_earth*0.15,...
        sin(nu_SOI_earth*0.85)*SOI_earth*0.15,...
        SOI_earth*0.07,0,'Color',col_sun,'LineWidth',1.5,'MaxHeadSize',0.7);
text(x_hyp_e(end),y_hyp_e(end)*1.1,z_hyp_e(end)+SOI_earth*0.09,...
     sprintf('v_{\\infty}=%.2f km/s',v_inf_earth),...
     'Color',col_sun,'FontSize',7,'FontName','Consolas');

view(-42,20);
xlabel('km','FontName','Consolas'); ylabel('km','FontName','Consolas');
zlabel('km','FontName','Consolas');
title('Phase 1 — Earth departure','Color','w','FontSize',10,'FontName','Consolas');

%% ============================================================
%  BOTTOM RIGHT: Mars arrival hyperbola
%% ============================================================
ax3 = subplot(2,3,6);
set(ax3,'Color',ax_bg,...
    'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5],'ZColor',[0.5 0.5 0.5],...
    'GridColor',[0.15 0.15 0.22],'GridAlpha',1,...
    'XGrid','on','YGrid','on','ZGrid','on',...
    'FontSize',8,'FontName','Consolas');
hold on; axis equal;

[xm3,ym3,zm3] = sphere_fn(R_mars, 0, 0, 0);
surf(xm3,ym3,zm3,'FaceColor',col_m,'EdgeColor','none','FaceAlpha',0.85);

plot3(r_park_mars*cos(theta_circ), r_park_mars*sin(theta_circ),...
      zeros(size(theta_circ)),'--','Color',[col_m 0.6],'LineWidth',1);
plot3(SOI_mars*cos(theta_circ), SOI_mars*sin(theta_circ),...
      zeros(size(theta_circ)),':','Color',[0.5 0.5 0.5 0.6],'LineWidth',0.7);
text(SOI_mars*0.62,SOI_mars*0.62,SOI_mars*0.05,'SOI',...
     'Color',[0.5 0.5 0.5],'FontSize',7,'FontName','Consolas');

z_hyp_m = linspace(SOI_mars*0.18, 0, length(x_hyp_m));
plot3(x_hyp_m, y_hyp_m, z_hyp_m,'Color',col_hy,'LineWidth',2);

scatter3(r_park_mars,0,0,80,col_bv,'filled','MarkerEdgeColor','w');
text(r_park_mars*1.1,r_park_mars*0.5,r_park_mars*0.2,...
     sprintf('MOI  \\DeltaV=%.2f',dv_MOI),...
     'Color',col_bv,'FontSize',7,'FontName','Consolas');

quiver3(x_hyp_m(1),y_hyp_m(1),z_hyp_m(1),...
        -cos(nu_SOI_mars*0.85)*SOI_mars*0.15,...
        -sin(nu_SOI_mars*0.85)*SOI_mars*0.15,...
        -SOI_mars*0.07,0,'Color',col_sun,'LineWidth',1.5,'MaxHeadSize',0.7);
text(x_hyp_m(1)*0.72,y_hyp_m(1)*1.1,z_hyp_m(1)+SOI_mars*0.07,...
     sprintf('v_{\\infty}=%.2f km/s',v_inf_mars),...
     'Color',col_sun,'FontSize',7,'FontName','Consolas');

view(-42,20);
xlabel('km','FontName','Consolas'); ylabel('km','FontName','Consolas');
zlabel('km','FontName','Consolas');
title('Phase 3 — Mars arrival','Color','w','FontSize',10,'FontName','Consolas');

sgtitle(sprintf('Patched Conic  |  Earth \\rightarrow Mars  |  \\DeltaV_{total} = %.3f km/s  |  TOF = %.0f days',...
        total_dv, T_days),...
        'Color','w','FontSize',13,'FontWeight','bold','FontName','Consolas');