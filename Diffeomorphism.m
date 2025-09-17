clc; clear; close all;

addpath(genpath('../'))

% W = 50; H = 50;
% 
% % Occupancy grid base
% M = -1*ones(H,W,'int8');   % -1 = unknown
% M(:,:) = 0;                % 0 = free
% 
% % Ostacolo 1: quadrato
% M(10:20, 10:20) = 100;
% 
% % Ostacolo 2: rettangolo
% M(30:40, 25:35) = 100;
% 
LAMBDA = 5e95; % LAMBDA should be big enough to ensure proper mapping


A = imread("map1.pgm");


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
for i=1:CF.NumObjects 

    current_b = b_reduced_history{i};
    if isempty(current_b) || size(current_b,1) < 3
        continue
    end

    b = double(current_b);
    b = b(all(isfinite(b),2),:);
    if size(b,1) < 3, continue; end

    % remove consecutive duplicate
    b = dedup_rc(b, 1e-9);
    if size(b,1) < 3, continue; end

    V   = [b(:,2), b(:,1)];

    polygon{i} = polyshape(V(:,1),V(:,2), 'Simplify', true, 'KeepCollinearPoints', true);
    polygon{i} = translate(polygon{i}, translation);
    polygon{i} = scale(polygon{i}, 0.05, [0 0]);

    current_b = polygon{i}.Vertices;
    current_b(size(current_b,1)+1,:) = current_b(1,:);

    vert{i} = current_b;

    [cx_p,cy_p] = centroid(polygon{i});
    ball_obstacle_centers{i} = [cx_p;cy_p];
    ball_obstacle_radius{i} = sqrt(area(polygon{i})/pi);
    real_obstacle_contours{i} = current_b.';

    realWorld.obstacles{end+1}.type    = 'qc';      
    realWorld.obstacles{end}.contour   = real_obstacle_contours{i};

    ballWorld.obstacles{end+1}.center        = ball_obstacle_centers{i};
    ballWorld.obstacles{end}.centerOriginal  = ball_obstacle_centers{i};
    ballWorld.obstacles{end}.radius          = ball_obstacle_radius{i};
    ballWorld.obstacles{end}.radiusOriginal  = ball_obstacle_radius{i};

    
end

figure;
for i=1:CF.NumObjects 
    % Set axis limits to match the original image size
    axis([realWorld.domain.contour(1,1) realWorld.domain.contour(1,4) realWorld.domain.contour(2,1) realWorld.domain.contour(2,2)]);
    axis equal;
    grid on;
    xlabel('X-pixels');
    ylabel('Y-pixels');
    set(gca,'YDir','reverse');
    title('Polygons');
    
    hold on;
    plot(polygon{i});
end
hold off;

% ===== Plot "mondo palla" (domain + ostacoli) =====
figure('Name','Ball world'); hold on; axis equal;
set(gca,'YDir','reverse');    % coord. immagine: y verso il basso

% Dominio (cerchio)
cB = ballWorld.domain.center(:);
RB = ballWorld.domain.radius;
th = linspace(0,2*pi,400);
plot(cB(1) + RB*cos(th), cB(2) + RB*sin(th), 'k','LineWidth',2);

% Ostacoli (cerchi)
Nobst = numel(ballWorld.obstacles);
for i = 1:Nobst
    ci = ballWorld.obstacles{i};
    Ri = ci.radius;
    plot(ci.center(1) + Ri*cos(th), ci.center(2) + Ri*sin(th), 'r','LineWidth',1.5);
    plot(ci.center(1), ci.center(2), 'r.','MarkerSize',14); % centro ostacolo
end

% Centro dominio e (opzionale) goal se coerente con le unità
plot(cB(1), cB(2), 'bo','MarkerFaceColor','b');         % centro ball
% Se il tuo goal è in stesse unità del ball, decommenta:
% gB = ballWorld.domain.goal(:);
% plot(gB(1), gB(2), 'gx','MarkerSize',10,'LineWidth',2);

% Cornice visiva
pad = 20;
xlim([cB(1)-RB-pad, cB(1)+RB+pad]);
ylim([cB(2)-RB-pad, cB(2)+RB+pad]);
grid on; title('Ball world (domain + obstacles)');


wm = WorldMapping(realWorld, ballWorld);
wm.evaluateMappings(LAMBDA);
[r2bMap, b2rMap, r2bJac, b2rJac] = wm.getMappings();

% --- costruisco f_k e beta_k come fa WorldMapping (solo per debug) ---
N = numel(realWorld.obstacles);
fH = cell(1,N+1); betaH = cell(1,N+1);

% ambiente (interior)
fH{1}    = ObstacleMappingQC(realWorld.domain.contour, 'interior', 0).getReal2BallMapHandle();
betaH{1} = @(q) 1 - norm(fH{1}(q))^2;

% ostacoli (exterior)
for k = 1:N
    fH{k+1}    = ObstacleMappingQC(realWorld.obstacles{k}.contour, 'exterior', 0).getReal2BallMapHandle();
    betaH{k+1} = @(q) norm(fH{k+1}(q))^2 - 1;
end

% prodotto delle beta escluso j
prodExcl = @(q,j) prod(arrayfun(@(k) betaH{k}(q), setdiff(1:N+1, j+1)));

% scansiona ogni ostacolo su un anello esterno sottile
lambda = LAMBDA;
viol = [];
for j = 1:N
    % anello: bordo ostacolo "just outside"
    ring = polyshape(realWorld.obstacles{j}.contour(1,:)', realWorld.obstacles{j}.contour(2,:)');
    % mezzo pixel fuori; regola se serve
    if isempty(ring), continue; end
    V = ring.Vertices; if isempty(V), continue; end
    V = V(1:10:end,:);                % campiona pochi punti

    minden = inf; minq = [NaN;NaN];
    for t = 1:size(V,1)
        q = V(t,:)';                   % [x;y]
        A = norm(q - realWorld.domain.goal)^2 * prodExcl(q, j);
        den = A + lambda * betaH{j+1}(q);
        if ~isfinite(den) || abs(den) < minden
            minden = abs(den); minq = q;
        end
    end

    if ~isfinite(minden) || minden < 1e-10
        fprintf('DEN≈0 per ostacolo %d  minden=%.3e  @ (%.3f, %.3f)\n', j, minden, minq(1), minq(2));
        viol(end+1) = j;
    end
end

if isempty(viol)
    disp('Nessun denominatore critico trovato sul bordo esterno degli ostacoli.');
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
