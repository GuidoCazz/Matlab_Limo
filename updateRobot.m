% Sketch of the Robot Pose update
function updateRobot(h, x, y, theta)
    theta = theta - pi/2;
    size = 0.2;
    pts = [ size, 0, -size;
            -size, 2*size, -size];
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
    ptsR = R*pts + [x;y];
    h.XData = ptsR(1,:);
    h.YData = ptsR(2,:);
end
