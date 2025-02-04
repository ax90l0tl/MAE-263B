clc
close all
clear
L1 = Link('revolute', 'd', 0.370, 'a', 0.0, 'alpha', 0.0, 'modified', 'qlim', [-deg2rad(170), deg2rad(170)]);
L2 = Link('revolute', 'd', 0.0, 'a', 0.325, 'alpha', 0.0, 'modified', 'qlim', [-deg2rad(145), deg2rad(145)]);
L3 = Link('revolute', 'd', 0.0, 'a', 0.225, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
L4 = Link('prismatic', 'a', 0, 'alpha', pi, 'modified', 'qlim', [0, 0.150]);
bot = SerialLink([L1, L2, L3, L4]);
% bot.teach();
% joint space
init_pose = [0.0, pi/2, 0.0, 0.0];
% world space
feeder_pos = [0.5, 0, 0.416 - 0.15/2, 0];
pcb_pos = [0.4, 0, 0.416 - 0.15/2, pi/4];
pcb_dim = [200, 100, 0]/1000;
chip_dim = [10, 10, 2]/1000; % l x w x h

pos_1 = [pcb_dim(1:2)/2 - chip_dim(1:2)/2, chip_dim(3), 0];
pos_2 = [-pcb_dim(1:2)/2 + chip_dim(1:2)/2, chip_dim(3), pi/2];
pos_3 = [pcb_dim(1)/2 - chip_dim(1)/2, -pcb_dim(2)/2 + chip_dim(2)/2, chip_dim(3), pi];
pos_4 = [-pcb_dim(1)/2 + chip_dim(1)/2, pcb_dim(2)/2 - chip_dim(2)/2, chip_dim(3), 3*pi/2];

pcb_pos_mat = transform(pcb_pos);
pos_1_mat = pcb_pos_mat * transform(pos_1);
pos_2_mat = pcb_pos_mat * transform(pos_2);
pos_3_mat = pcb_pos_mat * transform(pos_3);
pos_4_mat = pcb_pos_mat * transform(pos_4);
feeder_pos_mat = transform(feeder_pos);

pt = [0.2, 0.2, 0, 0.416];
% fk = IK(pcb_pos_mat);
% bot.fkine(fk)
% disp(fk)
% bot.plot(fk, 'workspace', [-1, 1, -1, 1, 0, 1], 'view', [0, 90]);

%%
clc
close all

a = 50*ones(1, 4) % deg/s^2
v1 = 420;
v2 = 720;
v3 = 1100;
v4 = 3000;
v = [v1, v2, v3, v4]; 

% fk = IK(pcb_pos);
fk = init_pose;
fk = [fk; IK(feeder_pos_mat, init_pose)];

linearParabolicBlendTrajectory2(fk, v, a, 10)

% fk = [fk; IK(pos_1_mat, fk(end, :))];
% fk = [fk; IK(feeder_pos_mat)];
% fk = [fk; IK(pos_2_mat)];
% fk = [fk; IK(feeder_pos_mat)];
% fk = [fk; IK(pos_3_mat)];
% fk = [fk; IK(feeder_pos_mat)];
% fk = [fk; IK(pos_4_mat)];

% disp(size(fk))
% fk
% jointspace_animation(fk, bot, length(fk), [0, 90], 'test')


%%

function joint_vals = IK(pt, prev_pose)
d4 = pt(3, 4);
a2 = 0.325;
a3 = 0.225;
l4 = 0.416 - 0.15;
% verify workspace
dist = pt(1, 4)^2 + pt(2, 4)^2;
theta = atan2(pt(2, 4), pt(1, 4));
if sqrt(dist) < a2 + a3 && sqrt(dist) > a2 - a3 && d4 <= 0.416 && d4 >= l4
    c_t2 = (a2^2 + a3^2 - dist)/(2*a2*a3);
    s_t2 = sqrt(1-c_t2^2);
    t2 = atan2(s_t2, c_t2) + pi;
    t2_alt = atan2(-s_t2, c_t2) + pi;
    s_t1 = (a3/sqrt(dist))*s_t2;
    c_t1 = sqrt(1-s_t1^2);
    t1 = atan2(s_t1, c_t1) + theta;
    t1_alt = atan2(-s_t1, c_t1) + theta;
    t3 = atan2(pt(2, 1), pt(1, 1)) - t2 - t1;
    t3_alt = -t3;
    % norm(abs([t1, t2, t3, d4] - prev_pose))
    % norm(abs([t1_alt, t2_alt, t3_alt, d4] - prev_pose))
    if norm(abs([t1, t2, t3, d4] - prev_pose)) > norm(abs([t1_alt, t2_alt, t3_alt, d4] - prev_pose))
        joint_vals = [t1_alt, t2_alt, t3_alt, d4];
    else
        joint_vals = [t1, t2, t3, d4];
    end
else
    disp('Location outside workspace')
    joint_vals = [NaN, NaN, NaN, NaN];
end
end

function T = transform(pt)
T = [[cos(pt(4)), -sin(pt(4)), 0, pt(1)];
    [sin(pt(4)), cos(pt(4)), 0, pt(2)];
    [0, 0, 1, pt(3)];
    [0, 0, 0, 1];];
end