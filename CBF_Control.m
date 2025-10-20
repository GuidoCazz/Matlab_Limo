clc; clear; close all

restoredefaultpath;
rehash toolboxcache;

qcm_path = fullfile(pwd, 'qcm-cbf');
ad_path = fullfile(pwd, 'MatlabAutoDiff');

addpath(genpath(qcm_path));
addpath(genpath(ad_path));

rng default;
rosIP   = '10.229.22.25';                 
MASTER  = 'http://10.229.22.87:11311';
CONTROL_OUT = 'cmd_vel';


LAMBDA = 5e10;
DT = 0.001;
T_MAX = 1e4;
POSITION_RADIUS = 'equal';
ALL_IN_BALL = false;
MAP_IN = 'map_raw';

rosshutdown; pause(0.2);
rosinit(MASTER,'NodeHost',rosIP,'NodeName','/matlab_control_final1');

tftree = rostf;
base_frame = 'base_link';    


matfile_path = 'diff_worldmapping.mat';
last_modified_time = 0;

pub = rospublisher(CONTROL_OUT,'geometry_msgs/Twist');
% Initialize the ROS node and set up the mapping process
disp('ROS Control:');

    
xTrajBall = zeros(2,T_MAX);
xTraj = zeros(2,T_MAX);

count = 0;
for t = 1 : T_MAX
    % Controlla se il matfile è stato modificato
    file_info = dir(matfile_path);
    if file_info.datenum > last_modified_time==0 || t==1

        % Se il file è stato modificato, aggiorna le mappature
        m = matfile(matfile_path);
        WM = m.wm;
        realWorld = m.realWorld;
        ballWorld = m.ballWorld;
    
        WM.evaluateMappings(LAMBDA);
        [r2bMap, b2rMap, r2bJac, b2rJac] = WM.getMappings();
    
        last_modified_time = file_info.datenum; % Aggiorna la data di modifica
        Nobst = numel(realWorld.obstacles);
        

        try
            tfBase = getTransform(tftree, 'map', 'base_link');            
            % Extract position
            x = [tfBase.Transform.Translation.X; tfBase.Transform.Translation.Y];

            q = [tfBase.Transform.Rotation.W, tfBase.Transform.Rotation.X, tfBase.Transform.Rotation.Y, tfBase.Transform.Rotation.Z];
            eul = quat2eul(q);       % default ZYX
            theta = eul(1);
        catch
            warning('TF not available at this timestep');
        end
        

        xBall = r2bMap(x);
        % xG = [4.5;-4.85];
        xG = [3;-2];
        xGBall = r2bMap(xG);
        

        xTraj(:,t) = x;
        xTrajBall(:,t) = xBall;

        % Plot iniziale
        figure;
        subplot(1,2,1), hold on, axis equal, set(gca, 'Visible', 'off')
        hSTraj = line(x(1), x(2), 'LineWidth', 2, LineJoin="chamfer");
        plot(realWorld.domain.contour(1,[1:end,1]), realWorld.domain.contour(2,[1:end,1]), 'LineWidth', 2, LineJoin="chamfer")
        for i = 1 : numel(realWorld.obstacles)
            plot(realWorld.obstacles{i}.contour(1,[1:end,1]), realWorld.obstacles{i}.contour(2,[1:end,1]), 'LineWidth', 2, LineJoin="chamfer")
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
            hObstBall{i} = plot(obstBall(1,[1:end,1]), obstBall(2,[1:end,1]), 'LineWidth', 2, LineJoin="chamfer");
        end
        hGBall = scatter(xGBall(1), xGBall(2), 200, '.');
        
        drawnow
    end
        
    % Continua con il ciclo di controllo
    if ALL_IN_BALL
        uNomBall = 10 * diag([6;1]) * (xGBall - xBall);
    else
        uNomReal = 10 * diag([6;1]) * (xG - x);
        Jr2b = AutoDiffJacobianFiniteDiff(r2bMap, x); % Calcola Jacobiano
        uNomBall = Jr2b * uNomReal; % Mappa in ball world
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
        Acbf(i,2*Nobst+i) = 2*ballWorld.obstacles{i}.radius;
        bcbf(i) = 2*(xBall-ballWorld.obstacles{i}.center)'*uBall + 0.5e2 * (norm(xBall-ballWorld.obstacles{i}.center)^2-ballWorld.obstacles{i}.radius^2);
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
        R = [cos(theta) -sin(theta);
             sin(theta)  cos(theta)];

%       u = r2bJac(x) \ uBall;
        u = Jr2b \ uBall;

        v = R' * u;
        v_max = 1.0;
        if abs(v(1)) > v_max
            v(1) = sign(v(1)) * v_max;
        end

        disp(v);
        d = [cos(theta); sin(theta)]; 
        dotp = dot(d,u);              
        crossp = d(1)*u(2) - d(2)*u(1);
        
        alpha = atan2(crossp, dotp);

        omega = 10*alpha;

        x = x + v(1)*[cos(theta);sin(theta)]*DT;
        theta = theta + omega*DT;
    end

    %%% plot
    for i = 1 : numel(ballWorld.obstacles)
        obstBall = ballWorld.obstacles{i}.center + ballWorld.obstacles{i}.radius * [cos(linspace(0,2*pi,100)); sin(linspace(0,2*pi,100))];
        hObstBall{i}.XData = obstBall(1,[1:end,1]);
        hObstBall{i}.YData = obstBall(2,[1:end,1]);
    end
    
    xTraj(:,t+1) = x;
    xTrajBall(:,t+1) = xBall;
    
    % hS.XData = x(1);
    % hS.YData = x(2);
    hSTraj.XData = xTraj(1,1:t+1);
    hSTraj.YData = xTraj(2,1:t+1);
    hSBall.XData = xBall(1);
    hSBall.YData = xBall(2);
    hSTrajBall.XData = xTrajBall(1,1:t+1);
    hSTrajBall.YData = xTrajBall(2,1:t+1);

    updateRobot(hRobot, x(1), x(2), theta);

    drawnow


    cmd = rosmessage(pub);

    cmd.Linear.X = v(1);
    cmd.Angular.Z = omega;

    send(pub, cmd);
end
