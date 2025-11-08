%% CAR LIKE VEHICLE ROS CONTROL
clc; clear; close all

restoredefaultpath;
rehash toolboxcache;

qcm_path = fullfile(pwd, 'qcm-cbf');
ad_path = fullfile(pwd, 'MatlabAutoDiff');

addpath(genpath(qcm_path));
addpath(genpath(ad_path));

rng default;
rosIP   = '129.97.71.90';                 
MASTER  = 'http://129.97.71.108:11311';
CONTROL_OUT = 'cmd_vel';

% Parameters Initialization
L = 0.2;
LAMBDA = 1e2; % Mapping Parameter
DT = 0.01;
T_MAX = 1e4;
POSITION_RADIUS = 'move';
ALL_IN_BALL = true;

rosshutdown; pause(0.2);
rosinit(MASTER,'NodeHost',rosIP,'NodeName','/matlab_control_final1');

tftree = rostf;
base_frame = 'base_link';    


matfile_path = 'diff_worldmapping.mat';
last_modified_time = 0;

pub = rospublisher(CONTROL_OUT,'geometry_msgs/Twist');
disp('ROS Control:');

%% First Matfile reading and Initialization
disp('Matfile arrived')
m = matfile(matfile_path);
WM = m.wm;
realWorld = m.realWorld;
ballWorld = m.ballWorld;

Nobst = numel(realWorld.obstacles);

% Lambda Dynamic Update
LAMBDA = 75*10^(Nobst);
if Nobst==1
    LAMBDA =75;
end

% Obstacle scaling factors
scale_width = 2.5; 
scale_length = 1; 

for i = 1 : Nobst
    contour = realWorld.obstacles{i}.contour;
    xi = contour(1, :);
    yi = contour(2, :);

    cx = mean(xi);
    cy = mean(yi);

    Xi = [xi - cx; yi - cy];

    C = cov(Xi');
    [V, ~] = eig(C);

    dir_long = V(:,2);
    dir_short = V(:,1);

    A = [dir_long dir_short] * diag([scale_length, scale_width]) * [dir_long dir_short]';

    X_scaled = A * Xi;

    realWorld.obstacles{i}.contour = [X_scaled(1,:) + cx; X_scaled(2,:) + cy];
end


WM.setRealWorld(realWorld)
WM.setBallWorld(ballWorld)
WM.evaluateMappings(LAMBDA)
[r2bMap, b2rMap, r2bJac, b2rJac] = WM.getMappings();
        

xTrajBall = zeros(2,T_MAX);
xTraj = zeros(2,T_MAX);

%% --- MAIN LOOP ---
for t = 1 : T_MAX
% Check for TF transformatiion
    try
        tfBase = getTransform(tftree, 'map', 'base_link');            
        x = [tfBase.Transform.Translation.X; tfBase.Transform.Translation.Y];

        q = [tfBase.Transform.Rotation.W, tfBase.Transform.Rotation.X, tfBase.Transform.Rotation.Y, tfBase.Transform.Rotation.Z];
        eul = quat2eul(q);
        theta = eul(1);
        x,theta
    catch
        warning('TF not available at this timestep');
    end

    xBall = r2bMap(x);
    xG = [5;0];
    xGBall = r2bMap(xG);    

    xTraj(:,t) = x;
    xTrajBall(:,t) = xBall;

    file_info = dir(matfile_path);

    % Obstacle Mapping Dynamic update
    if file_info.datenum > last_modified_time || t==1

        disp('Matfile arrived')
        m = matfile(matfile_path);
        WM = m.wm;
        realWorld = m.realWorld;
        ballWorld = m.ballWorld;
    
        Nobst = numel(realWorld.obstacles);
        LAMBDA = 75*10^(Nobst);
        if Nobst==1
            LAMBDA =75;
        end

        for i = 1 : Nobst
            contour = realWorld.obstacles{i}.contour;
            xi = contour(1, :);
            yi = contour(2, :);
        
            cx = mean(xi);
            cy = mean(yi);
        
            Xi = [xi - cx; yi - cy];
        
            C = cov(Xi');
            [V, ~] = eig(C);
        
            dir_long = V(:,2);
            dir_short = V(:,1);
        
            A = [dir_long dir_short] * diag([scale_length, scale_width]) * [dir_long dir_short]';
        
            X_scaled = A * Xi;
        
            realWorld.obstacles{i}.contour = [X_scaled(1,:) + cx; X_scaled(2,:) + cy];
        end

        WM.evaluateMappings(LAMBDA);
        [r2bMap, b2rMap, r2bJac, b2rJac] = WM.getMappings();
    
        last_modified_time = file_info.datenum;
        

        xBall = r2bMap(x);
        xGBall = r2bMap(xG);
        

        xTraj(:,t) = x;
        xTrajBall(:,t) = xBall;

        % --- Initial Plot ---
        obstacle_colors = ['r', 'g', 'b', 'c', 'm'];

        close all;
        figure;
        subplot(1,2,1), hold on, axis equal, set(gca, 'Visible', 'off')
        hSTraj = line(x(1), x(2), 'LineWidth', 2, LineJoin="chamfer");
        plot(realWorld.domain.contour(1,[1:end,1]), realWorld.domain.contour(2,[1:end,1]), 'LineWidth', 2, LineJoin="chamfer")
        for i = 1 : numel(realWorld.obstacles)
            plot(realWorld.obstacles{i}.contour(1,[1:end,1]), realWorld.obstacles{i}.contour(2,[1:end,1]), 'LineWidth', 2, LineJoin="chamfer",Color=obstacle_colors(i))
        end
        scatter(xG(1), xG(2), 200, '.');
        hRobot = drawRobot(x(1), x(2), theta, 0.2, 'r');


        subplot(1,2,2), hold on, axis equal, set(gca, 'Visible', 'off')
        hSBall = scatter(xBall(1), xBall(2), 200, '.');
        hSTrajBall = line(xBall(1), xBall(2), 'LineWidth', 2, LineJoin="chamfer");
        domainBall = ballWorld.domain.center + ballWorld.domain.radius * [cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))];
        plot(domainBall(1,[1:end,1]), domainBall(2,[1:end,1]), 'LineWidth', 2, LineJoin="chamfer")
        hObstBall = cell(1,numel(ballWorld.obstacles));
        for i = 1 : numel(ballWorld.obstacles)
            obstBall = ballWorld.obstacles{i}.center + ballWorld.obstacles{i}.radius * [cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))];
            hObstBall{i} = plot(obstBall(1,[1:end,1]), obstBall(2,[1:end,1]), 'LineWidth', 2, LineJoin="chamfer",Color=obstacle_colors(i));
        end
        hGBall = scatter(xGBall(1), xGBall(2), 200, '.');

        drawnow
    end
        
    % Nominal Controller
    if ALL_IN_BALL
        uNomBall = 10 * diag([1;1]) * (xGBall - xBall);
    else
        uNomReal = 10 * diag([1;1]) * (xG - x);
        Jr2b = AutoDiffJacobianFiniteDiff(r2bMap, x);
        uNomBall = Jr2b * uNomReal;
    end    


    %%% controller synthesis
    % robot does not care about the obstacles
    uBall = uNomBall;
    % obstacles care about the robot
    uObstNomBall = zeros(2*Nobst,1);
    rhoObstNomBall = zeros(Nobst,1);


    %% --- QP FORMULATION ---
    for i = 1 : Nobst
        uObstNomBall(2*i-1:2*i) = 10*(ballWorld.obstacles{i}.centerOriginal - ballWorld.obstacles{i}.center);
        rhoObstNomBall(i) = 10*(ballWorld.obstacles{i}.radiusOriginal - ballWorld.obstacles{i}.radius);
    end
    Acbf = zeros(Nobst,3*Nobst);
    bcbf = zeros(Nobst,1);
    for i = 1 : Nobst
        Acbf(i,2*i-1:2*i) = 2*(xBall-ballWorld.obstacles{i}.center)';
        Acbf(i,2*Nobst+i) = 2*(ballWorld.obstacles{i}.radius + L);
        bcbf(i) = 2*(xBall-ballWorld.obstacles{i}.center)'*uBall + 100 * (norm(xBall-ballWorld.obstacles{i}.center)^2-(ballWorld.obstacles{i}.radius + L)^2);
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

    %% Mappings re-evaluation
    WM.setBallWorld(ballWorld)
    WM.composeMappings(LAMBDA)
    [r2bMap, b2rMap, r2bJac, b2rJac] = WM.getMappings();

    %% Integration Step
    if ALL_IN_BALL
        xBall = xBall + uBall*DT;
        x_tr = b2rMap(xBall, x);

        d = [cos(theta); sin(theta)]; 
        direction = x_tr - x;

        v_proj = dot(d, direction) / (dot(d,d));
        v_max = 0.1;
        v = max(min(v_proj, v_max), -v_max);

        dotp = dot(d,direction);              
        crossp = d(1)*direction(2) - d(2)*direction(1);
        
        alpha = atan2(crossp, dotp);

        omega = 0.5*alpha;
    else
        Jr2b = AutoDiffJacobianFiniteDiff(r2bMap, x);%timeAD=toc;
        % R = [cos(theta) -sin(theta);
        %      sin(theta)  cos(theta)];
        d = [cos(theta); sin(theta)]; 

%       u = r2bJac(x) \ uBall;

        u = Jr2b \ uBall;

        v_proj = dot(d, u) / (dot(d,d))/500; 
        v_max = 0.1;
        v = max(min(v_proj, v_max), -v_max); 
        % if abs(uBall(1)) > v_max || abs(uBall(2)) > v_max
        %     uBall(1) = sign(uBall(1)) * v_max;
        %     uBall(2) = sign(uBall(2)) * v_max;
        % end

        dotp = dot(d,u);              
        crossp = d(1)*u(2) - d(2)*u(1);
        
        alpha = atan2(crossp, dotp);
        omega = 1.5*alpha;

        % if abs(v) < 0.05
        %     omega = 0.1*alpha;
        % else
        %     omega = 2*alpha;
        % end
        % x = x + v(1)*[cos(theta);sin(theta)]*DT;
        % theta = theta + omega*DT;

    end
    v,omega

    % Trajectories Plot 
    for i = 1 : numel(ballWorld.obstacles)
        obstBall = ballWorld.obstacles{i}.center + ballWorld.obstacles{i}.radius * [cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))];
        hObstBall{i}.XData = obstBall(1,[1:end,1]);
        hObstBall{i}.YData = obstBall(2,[1:end,1]);
    end
    
    xTraj(:,t+1) = x;
    xTrajBall(:,t+1) = xBall;
    
    hS.XData = x(1);
    hS.YData = x(2);
    hSTraj.XData = xTraj(1,1:t+1);
    hSTraj.YData = xTraj(2,1:t+1);
    hSBall.XData = xBall(1);
    hSBall.YData = xBall(2);
    hSTrajBall.XData = xTrajBall(1,1:t+1);
    hSTrajBall.YData = xTrajBall(2,1:t+1);

    updateRobot(hRobot, x(1), x(2), theta);
    drawnow limitrate

    %% Control Inputs Pubblication
    
    cmd = rosmessage(pub);

    cmd.Linear.X = v;
    cmd.Angular.Z = omega;

    send(pub, cmd);
end

