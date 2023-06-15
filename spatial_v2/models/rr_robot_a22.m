function  robot = rr_robot_a22
robot.NB = 2;
robot.parent = [0:robot.NB-1];

robot.ar{1} = [0;-1;0];
robot.ar{2} = [0;-1;0];

for i=1:robot.NB
    robot.L{i} = zeros(6,robot.NB);
    [~, S] = Rodrigues(robot.ar{i},0);
    robot.L{i}(:,i) = S;
end

% robot.Xtree{1} = rotx(pi/2);
% robot.Xtree{2} = xlt([0.0, -0.6, 0.0]);

robot.Xtree{1} = xlt([0.0, 0.0, 0.0]);
robot.Xtree{2} = xlt([0.0, 0.0, -1.0]);

% end effector
robot.kparent = [2];
robot.KXtree{1} = xlt([[0.0, 0.0, -0.7]]);

m1 = (0.3)*1; m2 = (0.3)*0.7;
I1 = (m1/12)*[(0.02)^2+(1)^2 0 0; 0 (0.05)^2+(1)^2 0; 0 0 (0.05)^2+(0.02)^2];
I2 = (m2/12)*[(0.02)^2+(0.7)^2 0 0; 0 (0.05)^2+(0.7)^2 0; 0 0 (0.05)^2+(0.02)^2];
robot.I{1} = mcI( m1, [0 0 -0.5], I1 ); 
robot.I{2} = mcI( m2, [0 0 -0.35], I2 ); 

robot.gravity = [0; 0; -9.81];

% oriented w.r.t Xtree{1}
% [x0 y0 z0; x1 y1 z1]
robot.appearance.body{1} = ...
    { 'box', [-0.07 0 -0.04; 0.07 -0.6 0.04], ...
      'cyl', [0 0 -0.07; 0 0 0.07], 0.1 };
  
robot.appearance.body{2} = ...
{ 'box', [-0.07 0 -0.04; 0.07 -0.4 0.04], ...
  'cyl', [0 0 -0.07; 0 0 0.07], 0.1 };