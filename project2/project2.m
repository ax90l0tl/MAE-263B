clc
close all
clear
L1 = Link('revolute', 'd', 0.416, 'a', 0.0, 'alpha', 0.0, 'modified', 'qlim', [-deg2rad(170), deg2rad(170)]);
L2 = Link('revolute', 'd', 0.0, 'a', 0.325, 'alpha', 0.0, 'modified', 'qlim', [-deg2rad(145), deg2rad(145)]);
L3 = Link('revolute', 'd', 0.0, 'a', 0.225, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
L4 = Link('prismatic', 'a', 0.0, 'alpha', 0.0, 'modified', 'qlim', [0, 0.150]);
bot = SerialLink([L1, L2, L3, L4]);


init_pose = [0.0, pi/2, 0.0, 0.0];
feeder_pos = [0.5, 0, 0.416 - 0.15/2, 0.0];
pcb_pos = [0.4, 0, 0.416 - 0.15/2, 0.0];


pt = [0.2, 0.2, 0, 0];
fk = IK(pcb_pos);
bot.fkine(fk)
disp(fk)


bot.plot(fk(2, :), 'workspace', [-1, 1, -1, 1, 0, 1], 'view', [0, 90]);

function joint_vals = IK(pt)
d4 = pt(4);
a2 = 0.325;
a3 = 0.225;
% verify workspace
dist = pt(1)^2 + pt(2)^2;
theta = atan2(pt(2), pt(1));
if sqrt(dist) < a2 + a3 && sqrt(dist) > a2 - a3
    c_t2 = (a2^2 + a3^2 - dist)/(2*a2*a3);
    s_t2 = sqrt(1-c_t2^2);
    t2 = atan2(s_t2, c_t2) + pi;
    t2_alt = atan2(-s_t2, c_t2) + pi;
    s_t1 = (a3/sqrt(dist))*s_t2;
    c_t1 = sqrt(1-s_t1^2);
    t1 = atan2(s_t1, c_t1) + theta;
    t1_alt = atan2(-s_t1, c_t1) + theta;
    t3 = pt(3) - t2 - t1;
    t3_alt = -t3;
else
    disp('Location outside workspace')
    t1 = 0;
    t2 = 0;
    t3 = 0;
end
joint_vals = [[t1, t2, t3, d4];
    [t1_alt, t2_alt, t3_alt, d4]];
end