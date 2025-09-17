clc; clear; close all;

addpath(genpath('../'))

rosIP   = '10.132.252.25';                 
MASTER  = 'http://10.132.252.87:11311';
MAP_IN  = '/map_raw';
MAP_OUT = '/map';

BYPASS  = false;            % true: copia 1:1; false: applica filtro
ACTIVATE_BAG = false;
BAG_READ = false;

ACKERMANN_RANGE = true;
BASE_FRAME = 'base_link';

% --- ROS ---
rosshutdown; pause(0.2);
rosinit(MASTER,'NodeHost',rosIP,'NodeName','/matlab_processing_map3');

tftree = rostf; pause(0.2);
pub = rospublisher(MAP_OUT,'nav_msgs/OccupancyGrid','IsLatching',true,'DataFormat','struct');
  

bagWriter = rosbagwriter('office_map_out.bag');
cleanup = onCleanup(@() close(bagWriter));

subProc = rossubscriber(MAP_IN,'nav_msgs/OccupancyGrid', @(~,m)processAndPublish(m,pub,BYPASS,bagWriter,ACTIVATE_BAG,BAG_READ,ACKERMANN_RANGE,tftree,BASE_FRAME), ...
    'DataFormat','struct');


disp('Relay attivo: /map_raw -> /map (latched).');

% ====== FUNZIONI ======
function processAndPublish(m,pub,BYPASS,bagWriter,ACTIVATE_BAG,BAG_READ,ACKERMANN_RANGE,tftree,BASE_FRAME)

    if ACTIVATE_BAG
        if BAG_READ
            bag = rosbag("office_map.bag");
            sel = select(bag,'Topic','/map_raw');
            msgs = readMessages(sel,'DataFormat','struct');
            if isempty(msgs), return; end
            m = msgs{end};
        else
            write(bagWriter,'/map_raw', m.Header.Stamp, m);    
        end
    end

    W = double(m.Info.Width);  H = double(m.Info.Height);
    if W==0 || H==0, return; end
    M = reshape(int16(m.Data), [W,H])';


    if ACKERMANN_RANGE
        % --- 0) Map info ---
        W   = double(m.Info.Width);
        H   = double(m.Info.Height);
        res = double(m.Info.Resolution);
        ox  = double(m.Info.Origin.Position.X);
        oy  = double(m.Info.Origin.Position.Y);
        map_frame = char(m.Header.FrameId);
    
        % --- 1) Vehicle params ---
        L          = 0.20;   % [m]
        Rmin       = 0.40;   % [m]
        kappa_max  = 1/Rmin; % [1/m]
        vmax       = 1.0;    % [m/s]
        Treach     = 3.0;    % [s]
        smax       = vmax * Treach;


        try
            % prendi la trasformazione più recente disponibile
            tr = getTransform(tftree, map_frame, BASE_FRAME);  % target, source
            lx = double(tr.Transform.Translation.X);
            ly = double(tr.Transform.Translation.Y);
            q  = tr.Transform.Rotation;
            yaw = atan2( 2*(q.W*q.Z + q.X*q.Y), 1 - 2*(q.Y*q.Y + q.Z*q.Z) );
        catch
            % fallback se TF non disponibile
            lx = ox; ly = oy; yaw = 0;
        end
        % --- 3) Frontiera (due lobi a curvatura satura) nel frame ROBOT ---
        Rmin      = 1/kappa_max;                 % raggio minimo
        theta_max = min(smax/Rmin, pi-1e-3);     % angolo massimo percorso sull'arco
        Ntheta    = 400;
        theta     = linspace(0, theta_max, Ntheta);
        
        % Archi sinistra/destra (frame robot, posa all'origine, heading +x)
        xL =  Rmin * sin(theta);
        yL =  Rmin * (1 - cos(theta));           % centro CL = (0, +Rmin)
        xR =  xL;
        yR = -yL;                                % centro CR = (0, -Rmin)
        
        % Settori (poligoni) chiusi: centro -> arco -> centro
        CL = [0, +Rmin];  CR = [0, -Rmin];
        PX_Lr = [CL(1),  xL, CL(1)];
        PY_Lr = [CL(2),  yL, CL(2)];
        PX_Rr = [CR(1),  xR, CR(1)];
        PY_Rr = [CR(2),  yR, CR(2)];
        
        % --- 4) Rototraslazione nel frame MAPPA con la posa corrente ---
        c = cos(yaw); s = sin(yaw);
        rot = @(xr,yr) deal(lx + c*xr - s*yr, ly + s*xr + c*yr);
        [PX_L, PY_L] = rot(PX_Lr, PY_Lr);
        [PX_R, PY_R] = rot(PX_Rr, PY_Rr);
        
        % --- 5) Maschera su griglia mappa ---
        [Xc,Yc] = meshgrid(1:W, 1:H);
        Xw = ox + (Xc - 0.5)*res;   % centri celle (m)
        Yw = oy + (Yc - 0.5)*res;
            
        % === maschera "luna" interna: due settori a raggio Rmin (già calcolati) ===
        inL = inpolygon(Xw, Yw, PX_L, PY_L);
        inR = inpolygon(Xw, Yw, PX_R, PY_R);
        inInner = inL | inR;                     % zona NON raggiungibile (limite sterzo)
        
        % === settore GRANDE in avanti (inReachable "lordo") ===
        Rbig = smax;                             % raggio max raggiungibile in 3 s
        phi  = min(pi/2, 1.2*(smax/Rmin) + 0.2); % apertura (regolabile). Tipico: ~90°
        Narc = 600;
        ang  = linspace(-phi, +phi, Narc);
        % poligono del settore nel frame robot
        PX_big_r = [0, Rbig*cos(ang), 0];
        PY_big_r = [0, Rbig*sin(ang), 0];
        % rototraslazione nel frame mappa
        [PX_big, PY_big] = rot(PX_big_r, PY_big_r);
        inBig = inpolygon(Xw, Yw, PX_big, PY_big);
        
        % === vincoli opzionali: solo avanti e piccolo "gap" vicino al robot ===
        dx = Xw - lx;  dy = Yw - ly;
        x_r =  c*dx + s*dy;    y_r = -s*dx + c*dy;
        inForward = (x_r >= 0);                % niente retro
        r_gap = 0.20;                          % buco vicino al robot (opzionale)
        notGap = (x_r.^2 + y_r.^2 >= r_gap^2);
        
        % === area raggiungibile = settore grande – luna interna, con vincoli ===
        inReachable = inBig & ~inInner & inForward & notGap;
        
        % applica alla mappa
        M(~inReachable) = -1;


    end


    
    if BYPASS
        M2 = M;
    else
        M2 = ProcessingDiffeomorphism(M);
    end

    out = rosmessage(pub);             
    out.Header = m.Header;
    t = rostime('now');
    out.Header.Stamp.Sec  = uint32(t.Sec);
    out.Header.Stamp.Nsec = uint32(t.Nsec);

    out.Info = m.Info;
    out.Data = int8(reshape(M2',[],1));


    % write(bagWriter,'/map', out.Header.Stamp, out);

    send(pub,out);
end

function M2 = ProcessingDiffeomorphism(M)
     
    A = zeros(size(M));

    A(M == 100) = 0;
    A(M == 0) = 254;
    A(M == -1) = 205;

    BW = imbinarize(A);
    BW_inverted = ~BW;
    
    BW_inverted = bwareaopen(BW_inverted, 3);
    se = strel('disk', 1);
    BW_inverted = imdilate(BW_inverted, se);
    
    CC = bwconncomp(BW_inverted, 8);
    p = regionprops(CC,"Area");
    areas = [p.Area];
    
    areaThreshold = 100;
    
    CH = false(size(BW_inverted));

    for i = 1:length(areas)
        pixelIdx = CC.PixelIdxList{i};
    
        tempBW = false(size(BW_inverted));
        tempBW(pixelIdx) = true;

        if areas(i) < areaThreshold
            tempCH = bwconvhull(tempBW, 'objects', 8);
            CH = CH | tempCH; 
        else
            CH = CH | tempBW;
        end
    end
    % ballWorld.domain.radius = sqrt(TotArea/pi);
    
    BW_filled = imfill(CH,"holes");
    boundaries = bwboundaries(CH);

    translation = [-1000,-1000];

    res = 0.05;
    
    CF = bwconncomp(CH, 8);
    tolerance = 0.03;


    b_reduced_history = cell(1, CF.NumObjects);
    for k=1:CF.NumObjects

       b = boundaries{k};
       b = reducepoly(b,tolerance);
       b = unique(b, 'rows','stable');

       b_reduced_history{k} = b;

    end

    realWorld.domain.type = 'qc';
    realWorld.domain.contour = [-12.5 -12.5 7.5 7.5 -12.5;-12.5 12.5 12.5 -12.5 -12.5];  
    
    ballWorld.domain.center = [0;0];
    realWorld.domain.goal = [0;2];
    
    ballWorld.domain.goal = realWorld.domain.goal;
    ballWorld.domain.radius = 12.50;

    realWorld.obstacles = {};
    ballWorld.obstacles = {};
  
    polygon = cell(1, CF.NumObjects);
    poly_plot = cell(1, CF.NumObjects);

    for i=1:CF.NumObjects 
    
        current_b = b_reduced_history{i};
        if isempty(current_b) || size(current_b,1) < 3, continue; end        
        
        b = double(current_b);
        b = b(all(isfinite(b),2),:);
        if size(b,1) < 3, continue; end
    
        % remove consecutive duplicate
        b = dedup_rc(b, 1e-9);
        if size(b,1) < 3, continue; end
    
        V   = [b(:,2), b(:,1)];
          
        polygon{i} = polyshape(V(:,1),V(:,2), 'Simplify', true, 'KeepCollinearPoints', true);
        poly_plot{i} = polygon{i};
        polygon{i} = translate(polygon{i}, translation);
        polygon{i} = scale(polygon{i}, res, [0 0]);

        current_b = polygon{i}.Vertices;
        current_b(size(current_b,1)+1,:) = current_b(1,:);
        
        [cx_p,cy_p] = centroid(polygon{i});
        ball_obstacle_centers{i} = [cx_p;cy_p];
        ball_obstacle_radius{i} = sqrt(area(polygon{i})/pi);
        real_obstacle_contours{i} = current_b.';
        
        obstacles{i} = struct('cx',cx_p,'cy',cy_p,'r',ball_obstacle_radius{i});

        realWorld.obstacles{end+1}.type    = 'qc';      
        realWorld.obstacles{end}.contour   = real_obstacle_contours{i};
    
        ballWorld.obstacles{end+1}.center        = ball_obstacle_centers{i};
        ballWorld.obstacles{end}.centerOriginal  = ball_obstacle_centers{i};
        ballWorld.obstacles{end}.radius          = ball_obstacle_radius{i};
        ballWorld.obstacles{end}.radiusOriginal  = ball_obstacle_radius{i};
    end
   
 
    wm = WorldMapping(realWorld, ballWorld);


    Diffeomorphism = matfile("diff_worldmapping.mat",'Writable',true);
    save('diff_worldmapping.mat','wm','realWorld','ballWorld');

    [height, width] = size(A);  
    pgmImage = false(height, width);
 

    for i = 1:length(polygon)
        if ~isempty(poly_plot{i})
        mask = createMask(poly_plot{i}, height, width);
        pgmImage = pgmImage | mask; 
        end
    end

    % ====== CONVERSIONE AL FORMATO OccupancyGrid ======

    M2 = -1*ones(size(A),'int8');

    M2(pgmImage ~= 0) = 100;
    M2(pgmImage == 0) = 0;

    % [X, Y] = meshgrid(1:width, 1:height);   % X: columns, Y: rows
    % 
    % % outside the ball domain => occupied
    % occ = (X - ballWorld.domain.center(1)).^2 + (Y - ballWorld.domain.center(2)).^2 > ((ballWorld.domain.radius)^2);
    % 
    % % obstacle discs => occupied 
    % for i = 1:numel(realWorld.obstacles)
    %     ci = obstacles{i};
    %     occ = occ | ((X - ci.cx).^2 + (Y - ci.cy).^2) <= (ci.r^2);
    % end
    % 
    % % assemble occupancy grid
    % M_ball = -1*ones(height, width, 'int8');   % unknown
    % M_ball(~occ) = 0;                 % free inside the ball and outside discs
    % M_ball( occ) = 100;               % occupied outside ball or inside discs

end


% Function to create a binary mask from a polyshape
function mask = createMask(polygon, height, width)
    [X, Y] = meshgrid(1:width, 1:height);
    mask = inpolygon(X, Y, polygon.Vertices(:,1), polygon.Vertices(:,2));
end


function b2 = dedup_rc(b, tol)
    % rimuove duplicati consecutivi e chiusura duplicata (input [row col])
    keep = true(size(b,1),1);
    for k = 2:size(b,1)
        if norm(b(k,:)-b(k-1,:)) <= tol
            keep(k) = false;
        end
    end
    b = b(keep,:);
    if size(b,1) >= 2 && norm(b(end,:)-b(1,:)) <= tol
        b(end,:) = [];
    end
    b2 = b;
end

