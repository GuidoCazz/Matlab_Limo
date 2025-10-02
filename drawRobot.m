function h = drawRobot(x, y, theta, size, color)
    % Triangolo isoscele con base size e altezza 2*size
    theta = theta - pi/2;
    pts = [ size, 0, -size;
            -size, 2*size, -size];
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
    ptsR = R*pts + [x;y];
    h = fill(ptsR(1,:), ptsR(2,:), color, 'FaceAlpha',0.3, 'EdgeColor',color, 'LineWidth',2);
end