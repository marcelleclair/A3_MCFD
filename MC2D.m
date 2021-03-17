set(0, 'DefaultFigureWindowStyle', 'docked');
close all

% Physical Constants
global C;
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

m_n = 0.26 * C.m_0;
T_i = 300; % initial temperature (K)
tau_mn = 0.2e-12; % given mean time between collisions
t = 0; % current time
n = 0; % current number of steps
dt = 5e-15; % time step duration
t_max = 1000*dt; % run for 1000 cycles
P_sca = 1 - exp(-dt / tau_mn); % scattering prob. per e- per dt
rho = 1e15; % electron densigy (cm^-2)
% P_sca = 0;

sig = sqrt((C.kb * T_i) / m_n); %std for maxwell-boltzmann
v_th = sqrt((2 * C.kb * T_i) / m_n); % mean thermal velocity of electrons
l_mn = v_th * tau_mn;

% Define spatial and temporal boundaries
num_e = 30000; % number of electrons
num_disp = 10; % number of electrons displayed
x_max = 200e-9; % maximum x position (nm)
y_max = 100e-9; % maximum y position (nm)

% Set static electric field
V0 = 0.1;
Ex = V0/x_max; % 500 kV/m due to 0.1 V drop over 200 nm
Ey = 0;
Fx = Ex*C.q_0;
Fy = Ey*C.q_0;

% boxes extend past outer bounds to top electrons from bouncing
% between the box and boundary
box1 = Obstruction([80e-9 -10e-9], 40e-9, 50e-9, 0);
box2 = Obstruction([80e-9 60e-9], 40e-9, 50e-9, 0);
boxes = [box1 box2]; % array of all obstructions

% Find potential with resistive boxes included
% V0 = 0.1;
% cMap = makeCMap(x_max,y_max,2e-9,boxes);
% [V,Ex,Ey] = v_fd(V0,cMap,2e-9);

% t_slc = zeros(1, num_e); % Times since last colision
% col_total = 0; % cumulative number of electron collisions and scatters
% MFP = 0;
% tau_calc = 0;
% t_free = []; % list of all times between collisions
% l_free = []; % list of all free path lengths

% randomly distribute electrons within the boundaries
x = x_max * rand(1, num_e);
xp = x;
y = y_max * rand(1, num_e);
yp = y;

% don't let electrons start in the boxes
in_box = 1;
while(any(in_box))
    in_box = false(1, num_e);
    for b = 1 : length(boxes)
        in_box = in_box | isInside([x;y], boxes(b));
    end
    x(in_box) = x_max * rand(1, nnz(in_box));
    xp(in_box) = x(in_box);
    y(in_box) = y_max * rand(1, nnz(in_box));
    yp(in_box) = y(in_box);
end
% assign random speed by maxwell-boltzmann
vx = normrnd(0,sig,1,num_e);
vy = normrnd(0,sig,1,num_e);

% % assign constant speed and random direction
% theta = 2 * pi * rand(1, num_e);
% vx = v_th * cos(theta);
% vy = v_th * sin(theta);

% temp as function of electron velocities
temp = @(vx, vy) (m_n .* mean((vx.^2) + (vy.^2))) / (2 * C.kb);
T = temp(vx, vy);
Tp = T; % temperature at last step
T_avg = T;

% current as function of electron velocities
Ix = @(vx) mean(vx)*C.q_0*rho*(10e4)*y_max;
Iy = @(vy) mean(vy)*C.q_0*rho*(10e4)*x_max;
Ix_v = [];
t_v = [];

% initialize figures
fig_traj = figure("Name", "Trajectories");
ax_traj = gca;
pbaspect([2 1 1]);
xlim(ax_traj, [0 x_max]);
ylim(ax_traj, [0 y_max]);
hold(ax_traj, 'on');
% draw boxes
for b = 1 : length(boxes)
    rectangle('Position', [boxes(b).origin boxes(b).x_size boxes(b).y_size]);
end
fig_temp = figure("Name", "Temperature");
ax_temp = gca;
title(ax_temp, "Semiconductor Temperature with Maxwell-Boltzmann Scattering");
ylabel(ax_temp, "Temperature (K)");
xlabel(ax_temp, "Time (s)");
hold(ax_temp, 'on');
grid(ax_temp, 'on');
fig_hist = figure("Name", "Velocity Distribution");
ax_hist = gca;
fig_I = figure("Name", "Current");
ax_I = gca;
title(ax_I, "Current (V_0 = " + V0 + ")");
xlabel(ax_I, "Time (s)");
ylabel(ax_I, "I_x (A)");
%hold(ax_I, 'on');
grid(ax_I, 'on');

while t < t_max
    % compute new position
%     x = x + vx * dt;
%     y = y + vy * dt;
    x = xp + vx*dt + 0.5*Fx*(dt^2)/m_n;
    vx = (x - xp) ./ dt;
    y = yp + vy*dt + 0.5*Fy*(dt^2)/m_n;
    vy = (y - yp) ./ dt;
    % scatter electrons
    scatter = rand(1,num_e) < P_sca;
    vx(scatter) = normrnd(0,sig,1,nnz(scatter));
    vy(scatter) = normrnd(0,sig,1,nnz(scatter));
    % detect collisions and apply boundary conditions
    x_collision = (x < 0) | (x > x_max);
    x(x_collision) = mod(x(x_collision), x_max);
    xp(x_collision) = x(x_collision);
%     y_collision = (y < 0) | (y > y_max);
%     vy(y_collision) = -vy(y_collision);
    yc_top = (y > y_max);
    vy(yc_top) = -abs(vy(yc_top));
    yc_bot = (y < 0);
    vy(yc_bot) = abs(vy(yc_bot));
    % handle collisions with obstacles
    for b = 1 : length(boxes)
        [IN, INL, INR, INU, IND] = isInside([x;y], boxes(b));
        if boxes(b).BC == 0
            % specular boundary
            vx(INL) = -abs(vx(INL));
            vx(INR) = abs(vx(INR));
            vy(IND) = -abs(vy(IND));
            vy(INU) = abs(vy(INU));
        elseif boxes(b).BC == 1
            vx(INL) = -abs(normrnd(0,sig,1,nnz(scatter)));
            vx(INR) = abs(normrnd(0,sig,1,nnz(scatter)));
            vy(IND) = -abs(normrnd(0,sig,1,nnz(scatter)));
            vy(INU) = abs(normrnd(0,sig,1,nnz(scatter)));
        end
    end
    % update temp and velocity magnitude, v is left to now
    % to allow computing mean-free-path
    T = temp(vx, vy);
    v = sqrt(vx.^2 + vy.^2);
    Ix_v = [Ix_v Ix(vx)];
    t_v = [t_v t];
    % advance clock
    t = t + dt;
    n = n + 1;
    T_avg = T_avg + (T - T_avg)/n; % rolling average temperature
    fprintf("Time: %3.3E s / %3.3E s\n", t, t_max);
    title(ax_traj, "Temperature = " + T + " K");
    set(ax_traj, 'ColorOrderIndex',1);
    plot(ax_traj, [xp(1:num_disp); x(1:num_disp)], [yp(1:num_disp); y(1:num_disp)]);
    plot(ax_temp, [(t-dt) t], [Tp T], 'r');
    plot(ax_temp, t, T_avg, '.g');
    histogram(ax_hist, v, 'Normalization', 'probability');
    title(ax_hist, "Electron Velocity Distribution");
    xlabel(ax_hist, "Velocity (m/s)");
    ylabel(ax_hist, "Relative Frequency");
    plot(ax_I, t_v, Ix_v);
    pause(0.00001); 
    % present becomes past
    Tp = T;
    xp = x;
    yp = y;
end

n_bins_y = 20;
n_bins_x = n_bins_y * 2;
A_bin = (100 * y_max / n_bins_y)^2; % area of each square in cm^2
edges_y = linspace(0, y_max, n_bins_y + 1);
edges_x = linspace(0, x_max, 2*n_bins_y + 1);

bin_y = discretize(y, edges_y);
bin_x = discretize(x, edges_x);

hmap = zeros(n_bins_y, n_bins_x);
dens = zeros(n_bins_y, n_bins_x);

for ni = 1 : n_bins_x
    for nj = 1 : n_bins_y
        hmap(nj, ni) = temp(vx(bin_x == ni & bin_y == nj), vy(bin_x == ni & bin_y == nj));
        dens(nj, ni) = nnz(bin_x == ni & bin_y == nj) .*((1e4*x_max*y_max*rho)/(num_e * A_bin));
    end
end

fig_dens = figure("Name", "Electron Density");
surf(dens);
zlabel("Electron density (cm^{-3})");
xlabel('x (m)');
xticks([0 10 20 30 40]);
xticklabels([0 0.25 0.5 0.75 1]*x_max);
ylabel('y (m)');
yticks([0 10 20]);
yticklabels([0 0.5 1]*y_max);
pbaspect([2 1 1]);
view([1 -2 2]);

fig_heat = figure("Name", "Temperature Map");
pcolor(hmap);
annotation('textarrow',[1,1],[0.5,0.5],'string','Temperature (K)', ...
      'HeadStyle','none','LineStyle','none','HorizontalAlignment', ...
      'center','TextRotation',90);

xlabel('x (m)');
xticks([0 10 20 30 40]);
xticklabels([0 0.25 0.5 0.75 1]*x_max);
ylabel('y (m)');
yticks([0 10 20]);
yticklabels([0 0.5 1]*y_max);
colorbar
pbaspect([2 1 1]);
 
% hist3(transpose([x;y]), [20,10]);
% pbaspect([2 1 1]);
% view([1 -1 1]);

% fprintf("SIMULATION END AFTER %d STEPS\n", n);
% fprintf("Mean Free Path = %3.3E m\n", MFP);
% fprintf("Mean time between collisions = %3.3E\n", tau_calc);
