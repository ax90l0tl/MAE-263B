clc
close all
clear
L1 = Link('revolute', 'd', 0.400, 'a', 0.0, 'alpha', 0.0, 'modified', 'qlim', [-deg2rad(170), deg2rad(170)]);
L2 = Link('revolute', 'd', 0.0, 'a', 0.325, 'alpha', 0.0, 'modified', 'qlim', [-deg2rad(145), deg2rad(145)]);
L3 = Link('revolute', 'd', 0.0, 'a', 0.225, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
L4 = Link('prismatic', 'a', 0, 'alpha', pi, 'modified', 'qlim', [0, 0.150]);
bot = SerialLink([L1, L2, L3, L4]);
tool = transl(0, 0, 30/1000)*trotx(pi);
min_height = 0.4 - 0.15 - 30/1000;
max_height = 0.4 - 30/1000;
bot.tool = tool;

% bot.teach()

% joint space
init_pose = [0.0, pi/2, 0.0, 0.0];
% world space
feeder_pos_v1 = [0, 0.5, max_height, 0];
feeder_pos = [0, 0.5, min_height, 0];
pcb_pos = [-0.0, 0.35, min_height, 0];
pcb_dim = [100, 100, 0]/1000;
chip_dim = [10, 10, 0]/1000; % l x w x h

pos_1 = [pcb_dim(1:2)/2 - chip_dim(1:2)/2, chip_dim(3), 0];
pos_3 = [-pcb_dim(1:2)/2 + chip_dim(1:2)/2, chip_dim(3), pi];
pos_2 = [pcb_dim(1)/2 - chip_dim(1)/2, -pcb_dim(2)/2 + chip_dim(2)/2, chip_dim(3), pi/2];
pos_4 = [-pcb_dim(1)/2 + chip_dim(1)/2, pcb_dim(2)/2 - chip_dim(2)/2, chip_dim(3), 3*pi/2];

pcb_pos_mat = transform(pcb_pos)
pos_1_mat = pcb_pos_mat * transform(pos_1);
pos_1_v1_mat = pos_1_mat * transl(0, 0, 0.150);
pos_2_mat = pcb_pos_mat * transform(pos_2);
pos_2_v1_mat = pos_1_v1_mat * transl(0.02, 0.02, 0);
pos_2_v2_mat = pos_2_mat * transl(-0.01, -0.01, 0.15);
pos_2_v3_mat = pos_2_mat * transl(0, 0, 0.150);
pos_3_mat = pcb_pos_mat * transform(pos_3);
pos_4_mat = pcb_pos_mat * transform(pos_4);
feeder_pos_mat = transform(feeder_pos);
feeder_pos_v1_mat = transform(feeder_pos_v1);



bot.plot(init_pose, 'workspace', [-0.2, 0.6, -0.2, 0.6, 0, 1], 'view', [0, 90]);
hold on
trplot(feeder_pos_mat, 'frame', 'Feeder', 'length', 0.05, 'color', 'red', 'text_opts', {'FontSize', 8});
trplot(pcb_pos_mat, 'frame', 'PCB', 'length', 0.05, 'color', 'blue', 'text_opts', {'FontSize', 8});
trplot(pos_1_mat, 'frame', 'Pos1', 'length', 0.05, 'color', 'black', 'text_opts', {'FontSize', 6});
trplot(pos_2_v1_mat, 'frame', 'Pos2_v1', 'length', 0.05, 'color', 'black', 'text_opts', {'FontSize', 6});
trplot(pos_2_v2_mat, 'frame', 'Pos2_v2', 'length', 0.05, 'color', 'black', 'text_opts', {'FontSize', 6});
trplot(pos_2_mat, 'frame', 'Pos2', 'length', 0.05, 'color', 'black', 'text_opts', {'FontSize', 6});
trplot(pos_3_mat, 'frame', 'Pos3', 'length', 0.05, 'color', 'black', 'text_opts', {'FontSize', 6});
trplot(pos_4_mat, 'frame', 'Pos4', 'length', 0.05, 'color', 'black', 'text_opts', {'FontSize', 6});

%% Joint space trajectory

clc
close all

a = 50*ones(1, 4)*pi/180; % rad/s^2
v1 = 420*pi/180;
v2 = 720*pi/180;
v3 = 1100; % mm/s
v4 = 3000*pi/180; 
v = [v1, v2, v3, v4]; 

% trajectory 1
fk = init_pose;
% bot.fkine(IK(feeder_pos_v1_mat, fk(end, :)))
% bot.plot(IK(feeder_pos_v1_mat, fk(end, :)))
fk = [fk; IK(feeder_pos_v1_mat, fk(end, :))];
fk = [fk; IK(feeder_pos_mat, fk(end, :))];
[q1, qd1, qdd1, qddd1, time1] = linearParabolicBlendTrajectory3(fk, a, 10);


% trajectory 2
fk = IK(feeder_pos_mat, fk(end, :));
fk = [fk; IK(feeder_pos_v1_mat, fk(end, :))];
fk = [fk; IK(pos_1_v1_mat, fk(end, :))];
fk = [fk; IK(pos_1_mat, fk(end, :))];
[q2, qd2, qdd2, qddd2, time2] = linearParabolicBlendTrajectory3(fk, a, 10);
% trajectory 3
fk = IK(pos_1_mat, fk(end, :));
fk = [fk; IK(pos_1_v1_mat, fk(end, :))];
fk = [fk; IK(feeder_pos_v1_mat, fk(end, :))];
fk = [fk; IK(feeder_pos_mat, fk(end, :))];
[q3, qd3, qdd3, qddd3, time3] = linearParabolicBlendTrajectory3(fk, a, 10);
% trajectory 4
fk = IK(feeder_pos_mat, fk(end, :));
fk = [fk; IK(feeder_pos_v1_mat, fk(end, :))];
fk = [fk; IK(pos_2_v1_mat, fk(end, :))];
fk = [fk; IK(pos_2_v2_mat, fk(end, :))];
fk = [fk; IK(pos_2_v3_mat, fk(end, :))];
fk = [fk; IK(pos_2_mat, fk(end, :))];
% bot.plot(IK(pos_2_mat, fk(end, :)))
% bot.fkine(IK(pos_2_mat, fk(end, :)))
[q4, qd4, qdd4, qddd4, time4] = linearParabolicBlendTrajectory3(fk, a, 10);
% trajectory 5
fk = IK(pos_2_mat, fk(end, :));
fk = [fk; IK(pos_2_v3_mat, fk(end, :))];
fk = [fk; IK(pos_2_v2_mat, fk(end, :))];
fk = [fk; IK(pos_2_v1_mat, fk(end, :))];
fk = [fk; IK(feeder_pos_v1_mat, fk(end, :))];
fk = [fk; IK(feeder_pos_mat, fk(end, :))];
[q5, qd5, qdd5, qddd5, time5] = linearParabolicBlendTrajectory3(fk, a, 10);



q = [q1;q2;q3;q4;q5];
qd = [qd1;qd2;qd3;qd4;qd5];
qdd = [qdd1;qdd2;qdd3;qdd4;qdd5];
qddd = [qddd1;qddd2;qddd3;qddd4;qddd5];
time = [time1,time2 + time1(end),time3 + time2(end) + time1(end),time4 + time3(end) + time2(end) + time1(end),time5 + time4(end) + time3(end) + time2(end) + time1(end)];
jointspace_animation(q, bot, length(q), [200, 20], 'all')
% jointspace_animation(q1, bot, length(q1), [0, 90], 'init_to_feeder')
% jointspace_animation(q2, bot, length(q2), [0, 90], 'feeder_to_pos1')
% jointspace_animation(q3, bot, length(q3), [0, 90], 'pos1_to_feeder')
% figure()
% hold on
% trplot(feeder_pos_mat, 'frame', 'Feeder', 'length', 0.05, 'color', 'red', 'text_opts', {'FontSize', 8});
% trplot(pcb_pos_mat, 'frame', 'PCB', 'length', 0.05, 'color', 'blue', 'text_opts', {'FontSize', 8});
% trplot(pos_1_mat, 'frame', 'Pos1', 'length', 0.05, 'color', 'black', 'text_opts', {'FontSize', 6});
% trplot(pos_2_mat, 'frame', 'Pos2', 'length', 0.05, 'color', 'black', 'text_opts', {'FontSize', 6});
% jointspace_animation(q4, bot, length(q4), [0, 90], 'feeder_to_pos2')
% jointspace_animation(q5, bot, length(q5), [0, 90], 'pos2_to_feeder')

%%
% traj 1
% Plot positions
figure();
hold on
plot(time, q)
title('Joint Positions Over Time')
xlabel('Time (s)')
ylabel('Position (rad or m)')
legend('J1', 'J2', 'J3', 'J4')
hold off

% Plot velocities
figure();
hold on
plot(time, qd)
title('Joint Velocities Over Time')
xlabel('Time (s)')
ylabel('Velocity (rad/s or m/s)')
legend('J1', 'J2', 'J3', 'J4')
hold off

% Plot accelerations
figure();
hold on
plot(time, qdd)
title('Joint Accelerations Over Time')
xlabel('Time (s)')
ylabel('Acceleration (rad/s^2 or m/s^2)')
legend('J1', 'J2', 'J3', 'J4')
hold off

% Plot jerks
figure();
hold on
plot(time, qddd)
title('Joint Jerks Over Time')
xlabel('Time (s)')
ylabel('Jerk (rad/s^3 or m/s^3)')
legend('J1', 'J2', 'J3', 'J4')
hold off

%% task space
clc
close all

a = 50*ones(1, 4)*pi/180; % rad/s^2
v1 = 420*pi/180;
v2 = 720*pi/180;
v3 = 1100; % mm/s
v4 = 3000*pi/180; 
v = [v1, v2, v3, v4]; 

% trajectory 1
fk1 = untransform(double(bot.fkine(init_pose)));
fk1 = [fk1; untransform(feeder_pos_v1_mat)];
fk1 = [fk1; untransform(feeder_pos_mat)];
[q1, qd1, qdd1, qddd1, time1] = linearParabolicBlendTrajectory3(fk1, a, 10);
fk1 = zeros(size(q1));
for i = 1:1:length(q1)
    transform(q1(i, :))
    if i == 1
        fk1(i, :) = IK(transform(q1(i, :)), init_pose);
    else
        fk1(i, :) = IK(transform(q1(i, :)), fk1(end, :));
    end
end

% trajectory 2
fk2 = untransform(feeder_pos_mat);
fk2 = [fk2; untransform(feeder_pos_v1_mat)];
fk2 = [fk2; untransform(pos_1_v1_mat)];
fk2 = [fk2; untransform(pos_1_mat)];
[q2, qd2, qdd2, qddd2, time2] = linearParabolicBlendTrajectory3(fk2, a, 10);
fk2 = zeros(size(q2));
for i = 1:1:length(q2)
    if i == 1
        fk2(i, :) = IK(transform(q2(i, :)), fk1(end, :));
    else
        fk2(i, :) = IK(transform(q2(i, :)), fk2(end, :));
    end
end

% trajectory 3
fk3 = untransform(pos_1_mat);
fk3 = [fk3; untransform(pos_1_v1_mat)];
fk3 = [fk3; untransform(feeder_pos_v1_mat)];
fk3 = [fk3; untransform(feeder_pos_mat)];
[q3, qd3, qdd3, qddd3, time3] = linearParabolicBlendTrajectory3(fk3, a, 10);
fk3 = zeros(size(q3));
for i = 1:1:length(q3)
    if i == 1
        fk3(i, :) = IK(transform(q3(i, :)), fk2(end, :));
    else
        fk3(i, :) = IK(transform(q3(i, :)), fk3(end, :));
    end
end

% trajectory 4
fk4 = untransform(feeder_pos_mat);
fk4 = [fk4; untransform(feeder_pos_v1_mat)];
fk4 = [fk4; untransform(pos_2_v1_mat)];
fk4 = [fk4; untransform(pos_2_v2_mat)];
fk4 = [fk4; untransform(pos_2_v3_mat)];
fk4 = [fk4; untransform(pos_2_mat)];
[q4, qd4, qdd4, qddd4, time4] = linearParabolicBlendTrajectory3(fk4, a, 10);
fk4 = zeros(size(q4));
for i = 1:1:length(q4)
    if i == 1
        fk4(i, :) = IK(transform(q4(i, :)), fk2(end, :));
    else
        fk4(i, :) = IK(transform(q4(i, :)), fk4(end, :));
    end
end

% trajectory 5
fk5 = untransform(pos_2_mat);
fk5 = [fk5; untransform(pos_2_v3_mat)];
fk5 = [fk5; untransform(pos_2_v2_mat)];
fk5 = [fk5; untransform(pos_2_v1_mat)];
fk5 = [fk5; untransform(feeder_pos_v1_mat)];
fk5 = [fk5; untransform(feeder_pos_mat)];
[q5, qd5, qdd5, qddd5, time5] = linearParabolicBlendTrajectory3(fk5, a, 10);
fk5 = zeros(size(q5));
for i = 1:1:length(q5)
    if i == 1
        fk5(i, :) = IK(transform(q5(i, :)), fk3(end, :));
    else
        fk5(i, :) = IK(transform(q5(i, :)), fk5(end, :));
    end
end

fk = [fk1;fk2;fk3;fk4;fk5];
q = [q1;q2;q3;q4;q5];
qd = [qd1;qd2;qd3;qd4;qd5];
qdd = [qdd1;qdd2;qdd3;qdd4;qdd5];
qddd = [qddd1;qddd2;qddd3;qddd4;qddd5];
time = [time1,time2 + time1(end),time3 + time2(end) + time1(end),time4 + time3(end) + time2(end) + time1(end),time5 + time4(end) + time3(end) + time2(end) + time1(end)];
jointspace_animation(fk, bot, length(fk), [200, 20], 'Trajectory_space_videos/all')
% jointspace_animation(fk1, bot, length(fk1), [0, 90], 'Trajectory_space_videos/init_to_feeder')
% jointspace_animation(fk2, bot, length(fk2), [0, 90], 'Trajectory_space_videos/feeder_to_pos1')
% jointspace_animation(fk3, bot, length(fk3), [0, 90], 'Trajectory_space_videos/pos1_to_feeder')
% jointspace_animation(fk4, bot, length(fk4), [0, 90], 'Trajectory_space_videos/feeder_to_pos2')
% jointspace_animation(fk5, bot, length(fk5), [0, 90], 'Trajectory_space_videos/pos2_to_feeder')

%%
% Plot positions
figure();
hold on
plot(time, fk)
title('Joint Positions Over Time')
xlabel('Time (s)')
ylabel('Position (rad or m)')
legend('J1', 'J2', 'J3', 'J4')
hold off

figure();
hold on
plot(time, q)
title('Tool Positions Over Time')
xlabel('Time (s)')
ylabel('Position (rad or m)')
legend('X', 'Y', 'Z', 'Yaw')
hold off

% Plot velocities
figure();
hold on
plot(time, qd)
title('Tool Velocities Over Time')
xlabel('Time (s)')
ylabel('Velocity (rad/s or m/s)')
legend('X', 'Y', 'Z', 'Yaw')
hold off

% Plot accelerations
figure();
hold on
plot(time, qdd)
title('Tool Accelerations Over Time')
xlabel('Time (s)')
ylabel('Acceleration (rad/s^2 or m/s^2)')
legend('X', 'Y', 'Z', 'Yaw')
hold off

% Plot jerks
figure();
hold on
plot(time, qddd)
title('Tool Jerks Over Time')
xlabel('Time (s)')
ylabel('Jerk (rad/s^3 or m/s^3)')
legend('X', 'Y', 'Z', 'Yaw')
hold off

%%

function joint_vals = IK(pt, prev_pose)
d4 = 0.4 - 0.03 - pt(3, 4);
a2 = 0.325;
a3 = 0.225;
l4 = 0.4 - 0.15 - 0.03;
% verify workspace
dist = pt(1, 4)^2 + pt(2, 4)^2;
theta = atan2(pt(2, 4), pt(1, 4));
if sqrt(dist) < a2 + a3 && sqrt(dist) > a2 - a3
    c_t2 = (a2^2 + a3^2 - dist)/(2*a2*a3);
    s_t2 = sqrt(1-c_t2^2);
    t2 = pi + atan2(s_t2, c_t2);
    t2_alt = pi + atan2(-s_t2, c_t2);
    s_t1 = (a3/sqrt(dist))*s_t2;
    c_t1 = sqrt(1-s_t1^2);
    t1 = atan2(s_t1, c_t1) + theta;
    t1_alt = atan2(-s_t1, c_t1) + theta;
    t3 = atan2(pt(2, 1), pt(1, 1)) - t2 - t1;
    t3_alt = atan2(pt(2, 1), pt(1, 1)) - t2_alt - t1_alt;
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

function pt = untransform(T)
    pt = [T(1, 4), T(2, 4), T(3, 4), atan2(T(2, 1), T(1, 1))];
end