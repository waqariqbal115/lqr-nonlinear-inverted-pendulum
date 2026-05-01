%% MATLAB CODE
%% Inverted Pendulum Project
% When Optimal Control Fails:
% LQR vs LQG vs Nonlinear Reality in Inverted Pendulum Stabilization
%
% This script builds directly on the class LQR pendulum problem, then extends it by:
%   1) testing the same LQR gain on the nonlinear plant,
%   2) adding torque saturation,
%   3) adding process disturbance and measurement noise, and
%   4) replacing full-state feedback with a Kalman filter (LQG).
%
% Requires Control System Toolbox for lqr and lqe.

clear; clc; close all;
rng(3);

%% Model parameters (upright-equilibrium linearization from homework)
a = 2.943;
b = 1.0;
A = [0 1; a 0];
B = [0; b];
C = [1 0];                % only theta is measured

%% LQR design (same gain reused in all tests)
Q = diag([30 5]);
R = 1;
K = lqr(A,B,Q,R);

%% Kalman filter / LQG observer design
W = diag([0.5 2.0]);      % process-noise intensity used for filter design
V = 0.01;                 % measurement-noise intensity used for filter design
L = lqe(A,eye(2),C,W,V);

%% Simulation settings
T  = 10.0;
dt = 0.002;
t  = 0:dt:T;
theta0_deg = 25;
x0 = [deg2rad(theta0_deg); 0];
umax = 2.5;               % realistic actuator limit used in nonlinear tests

%% Compare three cases
% 1) Linear plant + full-state LQR (ideal benchmark)
outLinear = simulate_case(t, x0, A, B, C, K, L, a, ...
    'linear_lqr', inf, 0.0, 0.0);

% 2) Nonlinear plant + same LQR gain + saturation
outNL = simulate_case(t, x0, A, B, C, K, L, a, ...
    'nonlinear_lqr', umax, 0.0, 0.0);

% 3) Nonlinear plant + LQG + saturation + disturbance + noisy measurement
outLQG = simulate_case(t, x0, A, B, C, K, L, a, ...
    'nonlinear_lqg', umax, 0.4, deg2rad(1.5));

%% Quantitative summary
settleLinear = settling_time(t, rad2deg(outLinear.x(1,:)), 1.0);
settleNL     = settling_time(t, rad2deg(outNL.x(1,:)),     1.0);
settleLQG    = settling_time(t, rad2deg(outLQG.x(1,:)),    2.0);
rmseEstDeg   = sqrt(mean(rad2deg(outLQG.x(1,:) - outLQG.xhat(1,:)).^2));

fprintf('\n=== Key results ===\n');
fprintf('K = [%.4f  %.4f]\n', K(1), K(2));
fprintf('L = [%.4f; %.4f]\n', L(1), L(2));
fprintf('Settling time (linear, +/-1 deg):    %.3f s\n', settleLinear);
fprintf('Settling time (nonlinear, +/-1 deg): %.3f s\n', settleNL);
fprintf('Settling time (LQG, +/-2 deg):       %.3f s\n', settleLQG);
fprintf('Estimation RMSE (theta):             %.3f deg\n', rmseEstDeg);

%% Figure 1: Response comparison
figure('Color','w','Position',[80 80 1100 650]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(t, rad2deg(outLinear.x(1,:)), 'LineWidth', 2.8); hold on;
plot(t, rad2deg(outNL.x(1,:)),     'LineWidth', 2.8);
plot(t, rad2deg(outLQG.x(1,:)),    'LineWidth', 2.4);
yline(0,'k-','LineWidth',0.8);
grid on;
ylabel('\theta (deg)');
title('Local optimality vs nonlinear reality');
legend('Linear + full-state LQR', 'Nonlinear plant + same LQR', ...
       'Nonlinear + LQG + noise/disturbance', 'Location','northeast');
text(0.02, 0.93, sprintf('Initial condition: %d deg, q(0)=0', theta0_deg), ...
    'Units','normalized', 'FontSize', 11);

nexttile;
plot(t, outLinear.u, 'LineWidth', 2.8); hold on;
plot(t, outNL.u,     'LineWidth', 2.8);
plot(t, outLQG.u,    'LineWidth', 2.4);
yline( umax,'--', 'LineWidth',1.4);
yline(-umax,'--', 'LineWidth',1.4);
grid on;
xlabel('Time (s)');
ylabel('u');
legend('Linear + full-state LQR', 'Nonlinear plant + same LQR', ...
       'Nonlinear + LQG + noise/disturbance', 'Torque limit', '', ...
       'Location','best');

%% Figure 2: Observer performance
figure('Color','w','Position',[100 100 1100 650]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(t, rad2deg(outLQG.x(1,:)),    'k', 'LineWidth', 2.8); hold on;
plot(t, rad2deg(outLQG.xhat(1,:)), 'LineWidth', 2.3);
plot(t(1:25:end), rad2deg(outLQG.y(1:25:end)), '.', 'MarkerSize', 9);
grid on;
ylabel('\theta (deg)');
title('Observer performance under noisy sensing');
legend('True angle', 'Estimated angle', 'Measured angle', 'Location','northeast');

nexttile;
errDeg = rad2deg(outLQG.x(1,:) - outLQG.xhat(1,:));
plot(t, errDeg, 'LineWidth', 2.5); hold on;
yline(0,'k-','LineWidth',0.8);
grid on;
xlabel('Time (s)');
ylabel('Estimation error (deg)');

%% Figure 3: Region of attraction (approximate) + phase trajectories
figure('Color','w','Position',[120 120 1200 620]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Sweep initial conditions
thGrid = linspace(-85,85,35);
qGrid  = linspace(-4,4,29);
S = zeros(length(qGrid), length(thGrid));
for i = 1:length(qGrid)
    for j = 1:length(thGrid)
        outTmp = simulate_case(0:0.01:6, [deg2rad(thGrid(j)); qGrid(i)], ...
            A, B, C, K, L, a, 'nonlinear_lqr', umax, 0.0, 0.0);
        thetaEnd = rad2deg(outTmp.x(1,end));
        qEnd     = outTmp.x(2,end);
        bounded  = max(abs(outTmp.x(1,:))) < pi;
        S(i,j)   = (abs(thetaEnd) < 3.0) && (abs(qEnd) < 0.3) && bounded;
    end
end

nexttile;
imagesc(thGrid, qGrid, S);
set(gca,'YDir','normal');
xlabel('Initial angle \theta_0 (deg)');
ylabel('Initial rate q_0 (rad/s)');
title('Approximate region of attraction');
colormap(gca,[0.75 0.00 0.15; 0.00 0.48 0.20]);
cb = colorbar; cb.Ticks = [0 1]; cb.TickLabels = {'fails','stabilizes'};
hold on; plot(theta0_deg,0,'kx','MarkerSize',12,'LineWidth',2.5);
text(theta0_deg+4,0.2,'nominal test','FontSize',11,'Color','k');

nexttile; hold on; grid on;
trialAngles = [15 35 55 70];
for kTrial = 1:length(trialAngles)
    outTmp = simulate_case(0:0.005:6, [deg2rad(trialAngles(kTrial)); 0], ...
        A, B, C, K, L, a, 'nonlinear_lqr', umax, 0.0, 0.0);
    plot(rad2deg(outTmp.x(1,:)), outTmp.x(2,:), 'LineWidth', 2.5);
end
xline(0,'k-','LineWidth',0.8); yline(0,'k-','LineWidth',0.8);
xlabel('\theta (deg)');
ylabel('q (rad/s)');
title('Phase trajectories reveal local behavior');
legend('15 deg','35 deg','55 deg','70 deg','Location','northeast');



%% ---------------- Local functions ----------------
function out = simulate_case(t, x0, A, B, C, K, L, a, mode, umax, disturbanceScale, measSigma)
    dt = t(2) - t(1);
    x    = x0(:);
    xhat = [0; 0];

    nx = length(t);
    X    = zeros(2,nx);
    Xhat = zeros(2,nx);
    U    = zeros(1,nx);
    Y    = zeros(1,nx);

    for k = 1:nx
        tk = t(k);

        switch mode
            case {'linear_lqr','nonlinear_lqr'}
                u = -K*x;
            case 'nonlinear_lqg'
                u = -K*xhat;
            otherwise
                error('Unknown mode');
        end

        u = min(max(u,-umax),umax);
        d = disturbanceScale * (0.18*sin(1.4*tk) + 0.07*sin(4.3*tk));

        switch mode
            case 'linear_lqr'
                xdot = A*x + B*u;
            otherwise
                xdot = [x(2);
                        a*sin(x(1)) + u + d];
        end

        x = x + dt*xdot;
        y = C*x + measSigma*randn;

        if strcmp(mode,'nonlinear_lqg')
            xhatdot = A*xhat + B*u + L*(y - C*xhat);
            xhat    = xhat + dt*xhatdot;
        else
            xhat = x;
        end

        X(:,k)    = x;
        Xhat(:,k) = xhat;
        U(k)      = u;
        Y(k)      = y;
    end

    out.x    = X;
    out.xhat = Xhat;
    out.u    = U;
    out.y    = Y;
end

function ts = settling_time(t, y, tol)
    idx = find(abs(y) > tol, 1, 'last');
    if isempty(idx)
        ts = 0;
    elseif idx == length(t)
        ts = t(end);
    else
        ts = t(idx+1);
    end
end
