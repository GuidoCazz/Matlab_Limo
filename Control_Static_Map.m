clc; close all; clear;

restoredefaultpath;
rehash toolboxcache;

qcm_path = fullfile(pwd, 'qcm-cbf');
ad_path = fullfile(pwd, 'MatlabAutoDiff');

addpath(genpath(qcm_path));
addpath(genpath(ad_path));


rng default;

L = 0.2;
LAMBDA = 5e7;
DT = 0.001;
T_MAX = 1e4;
POSITION_RADIUS = 'radius';
ALL_IN_BALL = false;


m = matfile('diff_worldmapping.mat');
WM = m.wm;
realWorld = m.realWorld;
ballWorld = m.ballWorld;
% 
% realWorld.obstacles = realWorld.obstacles(3);
% ballWorld.obstacles = ballWorld.obstacles(3);
% 
% WM.setRealWorld(realWorld);
% WM.setBallWorld(ballWorld);
% WM.composeMappings(LAMBDA);

WM.evaluateMappings(LAMBDA);
[r2bMap, b2rMap, r2bJac, b2rJac] = WM.getMappings();

x = [-2;0];
theta = 0;
xBall = r2bMap(x);

xG = [6;0.2];
xGBall = r2bMap(xG);

xTraj = zeros(2,T_MAX);
xTraj(:,1) = x;
xTrajBall = zeros(2,T_MAX);
xTrajBall(:,1) = xBall;

Nobst = numel(realWorld.obstacles);

% plots
figure;

subplot(1,2,1), hold on, axis equal, set(gca, 'Visible', 'off')
hSTraj = line(x(1), x(2), 'LineWidth', 2, LineJoin="chamfer");
plot(realWorld.domain.contour(1,[1:end,1]), realWorld.domain.contour(2,[1:end,1]), 'LineWidth', 2,LineJoin="chamfer")
for i = 1 : Nobst
    plot(realWorld.obstacles{i}.contour(1,[1:end,1]), realWorld.obstacles{i}.contour(2,[1:end,1]), 'LineWidth', 2,LineJoin="chamfer")
end
scatter(xG(1), xG(2), 200, '.');
hRobot = drawRobot(x(1), x(2), theta, 0.2, 'r');

subplot(1,2,2), hold on, axis equal, set(gca, 'Visible', 'off')
hSBall = scatter(xBall(1), xBall(2), 200, '.');
hSTrajBall = line(xBall(1), xBall(2), 'LineWidth', 2, LineJoin="chamfer");
domainBall = ballWorld.domain.center + ballWorld.domain.radius * [cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))];
plot(domainBall(1,[1:end,1]), domainBall(2,[1:end,1]), 'LineWidth', 2, LineJoin="chamfer")
hObstBall = cell(1,numel(ballWorld.obstacles));
for i = 1 : Nobst
    obstBall = ballWorld.obstacles{i}.center + ballWorld.obstacles{i}.radius * [cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))];
    hObstBall{i} = plot(obstBall(1,[1:end,1]), obstBall(2,[1:end,1]), 'LineWidth', 2, LineJoin="chamfer");
end
hGBall = scatter(xGBall(1), xGBall(2), 200, '.');

drawnow


for t = 1 : T_MAX

    xBall = r2bMap(x);
    if ALL_IN_BALL
        uNomBall = 10 * diag([6;1]) * (xGBall - xBall);
    else
        uNomReal = 10 * diag([6;1]) * (xG - x);
        Jr2b = AutoDiffJacobianFiniteDiff(r2bMap, x);%timeAD=toc;
        
        %uNomBall = r2bJac(x) * uNomReal;
        uNomBall = Jr2b * uNomReal;
    end    


    %%% controller synthesis
    % robot does not care about the obstacles
    uBall = uNomBall;
    % obstacles care about the robot
    uObstNomBall = zeros(2*Nobst,1);
    rhoObstNomBall = zeros(Nobst,1);

    for i = 1 : Nobst
        uObstNomBall(2*i-1:2*i) = 10*(ballWorld.obstacles{i}.centerOriginal - ballWorld.obstacles{i}.center);
        rhoObstNomBall(i) = 10*(ballWorld.obstacles{i}.radiusOriginal - ballWorld.obstacles{i}.radius);
    end
    Acbf = zeros(Nobst,3*Nobst);
    bcbf = zeros(Nobst,1);
    for i = 1 : Nobst
        Acbf(i,2*i-1:2*i) = 2*(xBall-ballWorld.obstacles{i}.center)';
        Acbf(i,2*Nobst+i) = 2*(ballWorld.obstacles{i}.radius + L);
        bcbf(i) = 2*(xBall-ballWorld.obstacles{i}.center)'*uBall + 50 * (norm(xBall-ballWorld.obstacles{i}.center)^2-(ballWorld.obstacles{i}.radius + L)^2);
    end
    if strcmp(POSITION_RADIUS, 'move')
        Wcenter = 1;
        Wradius = 1e3;
    elseif strcmp(POSITION_RADIUS, 'radius')
        Wcenter = 1e3;
        Wradius = 1;
    elseif strcmp(POSITION_RADIUS, 'equal')
        Wcenter = 1;
        Wradius = 1;
    end
    if Nobst == 0
        uObstBallrhoObstBall = [uObstNomBall; rhoObstNomBall];
    else
        uObstBallrhoObstBall = quadprog(2*blkdiag(Wcenter*eye(2*Nobst), Wradius*eye(Nobst)), ...
            -2*[Wcenter*uObstNomBall; Wradius*rhoObstNomBall]', ...
            Acbf, ...
            bcbf, ...
            [],[],[],[],[],optimoptions(@quadprog,'Display','off'));
    end
    for i = 1 : Nobst
        xDotI = uObstBallrhoObstBall(2*i-1:2*i);
        rDotI = uObstBallrhoObstBall(2*Nobst+i);
        ballWorld.obstacles{i}.center = ballWorld.obstacles{i}.center + xDotI*DT;
        ballWorld.obstacles{i}.radius = ballWorld.obstacles{i}.radius + rDotI*DT;
    end

    %%% re-evaluate mappings after changing obstacles
    WM.setBallWorld(ballWorld);
    WM.composeMappings(LAMBDA);
    [r2bMap, b2rMap, r2bJac, b2rJac] = WM.getMappings();

    %%% integration step
    if ALL_IN_BALL
        xBall = xBall + uBall*DT;
        x = b2rMap(xBall, x);
    else
        Jr2b = AutoDiffJacobianFiniteDiff(r2bMap, x);%timeAD=toc;
        % R = [cos(theta) -sin(theta);
        %      sin(theta)  cos(theta)];
        d = [cos(theta); sin(theta)]; 

%       u = r2bJac(x) \ uBall;

        u = Jr2b \ uBall;

        v_proj = dot(d, u) / (dot(d,d));   % scala (può essere negativa)
        v_max = 20;
        v = max(min(v_proj, v_max), -v_max);  % clamp preservando segno
        % if abs(uBall(1)) > v_max || abs(uBall(2)) > v_max
        %     uBall(1) = sign(uBall(1)) * v_max;
        %     uBall(2) = sign(uBall(2)) * v_max;
        % end

        dotp = dot(d,u);              
        crossp = d(1)*u(2) - d(2)*u(1);
        
        alpha = atan2(crossp, dotp);

        omega = 100*alpha;

        x = x + v*[cos(theta);sin(theta)]*DT;

        theta = theta + omega*DT;
        v,omega

    end

    %%% plot
    for i = 1 : numel(ballWorld.obstacles)
        obstBall = ballWorld.obstacles{i}.center + ballWorld.obstacles{i}.radius * [cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))];
        hObstBall{i}.XData = obstBall(1,[1:end,1]);
        hObstBall{i}.YData = obstBall(2,[1:end,1]);
    end
    
    xTraj(:,t) = x;
    xTrajBall(:,t) = xBall;

    
    % hS.XData = x(1);
    % hS.YData = x(2);
    hSTraj.XData = xTraj(1,1:t);
    hSTraj.YData = xTraj(2,1:t);
    hSBall.XData = xBall(1);
    hSBall.YData = xBall(2);
    hSTrajBall.XData = xTrajBall(1,1:t);
    hSTrajBall.YData = xTrajBall(2,1:t);

    updateRobot(hRobot, x(1), x(2), theta);

    drawnow
end

