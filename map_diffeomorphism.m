clc; clear; close all;

addpath(genpath('../'))

rosIP   = '10.71.141.25';                 
MASTER  = 'http://10.71.141.87:11311';
MAP_IN  = '/map_raw';
MAP_OUT = '/map';

BYPASS  = false;            % true: copia 1:1; false: applica filtro
ACTIVATE_BAG = false;
BAG_READ = false;

LAMBDA = 5e105;

% --- ROS ---
rosshutdown; pause(0.2);
rosinit(MASTER,'NodeHost',rosIP,'NodeName','/matlab_processing_map3');


pub = rospublisher(MAP_OUT,'nav_msgs/OccupancyGrid','IsLatching',true,'DataFormat','struct');
  

bagWriter = rosbagwriter('office_map_out.bag');
cleanup = onCleanup(@() close(bagWriter));

subProc = rossubscriber(MAP_IN,'nav_msgs/OccupancyGrid', @(~,m)processAndPublish(m,pub,BYPASS,bagWriter,ACTIVATE_BAG,BAG_READ,LAMBDA), ...
    'DataFormat','struct');

p_real = [1; 1];
if ~exist('R2BMAP','var')
    disp('R2BMAP non pronto: serve almeno un messaggio su /map_raw');
else
    p_ball = R2BMAP(p_real);
end



disp('Relay attivo: /map_raw -> /map (latched).');

% ====== FUNZIONI ======
function processAndPublish(m,pub,BYPASS,bagWriter,ACTIVATE_BAG,BAG_READ,LAMBDA)

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


    if BYPASS
        M2 = M; M_ball = M; r2bMap = [];
    else
        [M2, M_ball, r2bMap] = ProcessingDiffeomorphism(M,LAMBDA);
    end


    assignin('base','R2BMAP',r2bMap);

    out = rosmessage(pub);             
    out.Header = m.Header;
    t = rostime('now');
    out.Header.Stamp.Sec  = uint32(t.Sec);
    out.Header.Stamp.Nsec = uint32(t.Nsec);

    out.Info = m.Info;
    out.Data = int8(reshape(M_ball',[],1));


    % write(bagWriter,'/map', out.Header.Stamp, out);

    send(pub,out);
end

function [M2, M_ball, r2bMap] = ProcessingDiffeomorphism(M,LAMBDA)
     
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
    realWorld.domain.contour = [700 700 1200 1200 700;700 1200 1200 700 700];  
    
    ballWorld.domain.center = [1000;1000];
    realWorld.domain.goal = [1000;1040];
    
    ballWorld.domain.goal = realWorld.domain.goal;
    ballWorld.domain.radius = 300;

    realWorld.obstacles = {};
    ballWorld.obstacles = {};
  
    polygon = cell(1, CF.NumObjects);
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
        
        current_b = polygon{i}.Vertices;
        current_b(size(current_b,1)+1,:) = current_b(1,:);
        
        [cx_p,cy_p] = centroid(polygon{i});
        ball_obstacle_centers{i} = [cx_p;cy_p];
        ball_obstacle_radius{i} = sqrt(area(polygon{i})/pi);
        real_obstacle_contours{i} = current_b.';
        
        obstacles{i} = struct('cx',cx_p,'cy',cy_p,'r',ball_obstacle_radius);

        realWorld.obstacles{end+1}.type    = 'qc';      
        realWorld.obstacles{end}.contour   = real_obstacle_contours{i};
    
        ballWorld.obstacles{end+1}.center        = ball_obstacle_centers{i};
        ballWorld.obstacles{end}.centerOriginal  = ball_obstacle_centers{i};
        ballWorld.obstacles{end}.radius          = ball_obstacle_radius{i};
        ballWorld.obstacles{end}.radiusOriginal  = ball_obstacle_radius{i};

    end
 
    % wm.setRealWorld(realWorld); 
    % wm.setBallWorld(ballWorld);
    wm = WorldMapping(realWorld, ballWorld);
    wm.evaluateMappings(LAMBDA);
    [r2bMap, b2rMap, r2bJac, b2rJac] = wm.getMappings();

    [height, width] = size(A);  
    pgmImage = false(height, width);


    for i = 1:length(polygon)
        mask = createMask(polygon{i}, height, width);
        pgmImage = pgmImage | mask; 
    end

    % ====== CONVERSIONE AL FORMATO OccupancyGrid ======

    M2 = -1*ones(size(A),'int8');

    M2(pgmImage ~= 0) = 100;
    M2(pgmImage == 0) = 0;

    [X, Y] = meshgrid(1:width, 1:height);   % X: columns, Y: rows

    % outside the ball domain => occupied
    occ = (X - ballWorld.domain.center(1)).^2 + (Y - ballWorld.domain.center(2)).^2 > ((ballWorld.domain.radius)^2);

    % obstacle discs => occupied 
    for i = 1:Nobst
        ci = obstacles{i};
        occ = occ | ((X - ci.cx).^2 + (Y - ci.cy).^2) <= (ci.r^2);
    end

    % assemble occupancy grid
    M_ball = -1*ones(height, width, 'int8');   % unknown
    M_ball(~occ) = 0;                 % free inside the ball and outside discs
    M_ball( occ) = 100;               % occupied outside ball or inside discs
   
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
