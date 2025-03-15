clc
close all
clear

syms L1 L2 t1 t2 m1 m2

T_0_1 = trotx(0)*transl(0, 0, 0)*trotz(t1)*transl(0, 0, 0);
T_1_2 = trotx(0)*transl(L1, 0, 0)*trotz(t2)*transl(0, 0, 0);
T_2_tool = trotx(0)*transl(L2, 0, 0)*trotz(0)*transl(0, 0, 0);
T_0_2 = T_0_1*T_1_2;
T_0_tool = (T_0_1*T_1_2*T_2_tool);
R_0_1 = T_0_1(1:3, 1:3);
R_1_2 = T_1_2(1:3, 1:3);
R_2_tool = T_2_tool(1:3, 1:3);
R_0_2 = T_0_2(1:3, 1:3);
R_0_tool = T_0_tool(1:3, 1:3);
P_0_1 = T_0_1(1:3, end);
P_1_2 = T_1_2(1:3, end);
P_2_tool = T_2_tool(1:3, end);
P_0_2 = T_0_2(1:3, end);

%% Lagrangian
clc
syms I1xx I1yy I1zz I2xx I2yy I2zz
syms t1d t2d t1dd t2dd
syms g
I_1 = diag([I1xx, I1yy, I1zz]);
I_2 = diag([I2xx, I2yy, I2zz]);

I_1 = (1/12)*m1*L1^2*diag([0, 1, 1]);
I_2 = (1/12)*m2*L2^2*diag([0, 1, 1]);

I_0_1 = simplify(R_0_1*I_1*transpose(R_0_1));
I_0_2 = simplify(R_0_2*I_1*transpose(R_0_2));
P_c1 = [L1/2;0;0;1];
P_c2 = [L2/2;0;0;1];
P_0_c1 = T_0_1*P_c1;
P_0_c2 = simplify(T_0_2*P_c2);

Z_0_1 = T_0_1(1:3, 3);
Z_0_2 = T_0_2(1:3, 3);
J_0_v1 = [cross(Z_0_1, P_0_c1(1:3) - P_0_1), zeros(3, 1)]
J_0_w1 = [Z_0_1, zeros(3, 1)]

J_0_v2 = simplify([cross(Z_0_1, P_0_c2(1:3) - P_0_1), cross(Z_0_2, P_0_c2(1:3) - P_0_2)])
J_0_w2 = [Z_0_1, Z_0_2]

M = simplify(transpose(J_0_v1)*m1*J_0_v1 + transpose(J_0_w1)*I_0_1*J_0_w1 + transpose(J_0_v2)*m2*J_0_v2 + transpose(J_0_w2)*I_0_2*J_0_w2)

V_1 = sym(0);
V_2 = sym(0);
dq = [t1d, t2d];
th = [t1, t2];
for j = 1:2
    for k = 1:2
        V_1 = V_1 + (diff(M(1, j), th(k)) - 0.5*diff(M(j, k), t1))*dq(j)*dq(k);
        V_2 = V_2 + (diff(M(2, j), th(k)) - 0.5*diff(M(j, k), t2))*dq(j)*dq(k);
    end
end
V = simplify([V_1;V_2])

grav = [0; -g; 0];
G_1 = simplify(-m1*transpose(grav)*J_0_v1(:, 1) - m2*transpose(grav)*J_0_v2(:, 1));
G_2 = simplify(-m1*transpose(grav)*J_0_v1(:, 2) - m2*transpose(grav)*J_0_v2(:, 2));
G = [G_1;G_2]

Tau = simplify(M*[t1dd;t2dd] + V + G)
Tau(1)
Tau(2)

%% Newton Euler
clc
P_c1 = [L1/2;0;0];
P_c2 = [L2/2;0;0];

I_1 = diag([I1xx, I1yy, I1zz]);
I_2 = diag([I2xx, I2yy, I2zz]);

% I_1 = (1/12)*m1*L1^2*diag([0, 1, 1]);
% I_2 = (1/12)*m2*L2^2*diag([0, 1, 1]);

w_0 = zeros(3, 1);
wd_0 = zeros(3, 1);

v_0 = zeros(3, 1);
vd_0 = [0;g;0;];

f3 = [-10;0;0]; % N
n3 = [0;0;10]; % N*m
% Outward Iteration
%  i = 1
w_1 = transpose(R_0_1)*w_0 + t1d*[0;0;1];
wd_1 = transpose(R_0_1)*wd_0 + cross(transpose(R_0_1)*w_0, t1d*[0;0;1]) + t1dd*[0;0;1];
vd_1 = transpose(R_0_1)*(cross(wd_0, P_0_1) + cross(w_0, cross(w_0, P_0_1)) + vd_0);
vd_1_c = cross(wd_1, P_c1) + cross(w_1, cross(w_1, P_c1)) + vd_1;
F_1 = m1*vd_1_c;
N_1 = I_1*wd_1 + cross(w_1, I_1*w_1);


w_2 = transpose(R_1_2)*w_1 + t2d*[0;0;1];
wd_2 = transpose(R_1_2)*wd_1 + cross(transpose(R_1_2)*w_1, t2d*[0;0;1]) + t2dd*[0;0;1];
vd_2 = transpose(R_1_2)*(cross(wd_1, P_1_2) + cross(w_1, cross(w_1, P_1_2)) + vd_1);
vd_2_c = cross(wd_2, P_c2) + cross(w_2, cross(w_2, P_c1)) + vd_2;
F_2 = m2*vd_2_c;
N_2 = I_2*wd_2 + cross(w_2, I_2*w_2);

% Inward iteration
% i = 2
f2 = R_2_tool*f3 + F_2;
n2 = N_2 + R_2_tool*n3 + cross(P_c2, F_2) + cross(P_2_tool, R_2_tool*f3);
% i = 1
f1 = R_1_2*f2 + F_1;
n1 = N_1 + R_1_2*n2 + cross(P_c1, F_1) + cross(P_1_2, R_1_2*f2);

tau1 = simplify(transpose(n1)*[0, 0, 1]')
tau2 = simplify(transpose(n2)*[0, 0, 1]')
Tau = simplify([tau1;tau2]);
syms t1d2 t2d2
Tau = subs(Tau, [t1d^2, t2d^2], [t1d2, t2d2]);
%% c
clc
% close all
r1 = 0.1/2;
r2 = r1 - 0.005;
m_o = pi*(r1^2)*0.5*2710; % kg
m_i = pi*(r2^2)*0.5*2710;
mass = m_o - m_i; % kg
I = diag([0.5*m_o*r1^2, (1/12)*m_o*(3*r1^2 + 0.5^2), (1/12)*m_o*(3*r1^2 + 0.5^2)]) - diag([0.5*m_i*r2^2, (1/12)*m_i*(3*r2^2 + 0.5^2), (1/12)*m_i*(3*r2^2 + 0.5^2)]);
Tau = subs(Tau, [L1, L2, m1, m2], [0.5, 0.5, mass, mass]);
Tau = subs(Tau, [I1xx, I1yy, I1zz, I2xx, I2yy, I2zz], [I(1, 1), I(2, 2), I(3, 3), I(1, 1), I(2, 2), I(3, 3)]);

pos1 = ik(0.5, 0.5, [0.1, 0.0]);
pos2 = ik(0.5, 0.5, [0.9, 0.0]);

n = 20;
[q1, qd1, qdd1, qddd1, time1] = linearParabolicBlendTrajectory(0, 4, pos1(1), pos2(1), 0, 5, n);
[q2, qd2, qdd2, qddd2, time2] = linearParabolicBlendTrajectory(0, 4, pos1(2), pos2(2), 0, 5, n);
time = time1;
q = transpose([q1; q2]);
qd = transpose([qd1; qd2]);
qdd = transpose([qdd1; qdd2]);

L(1) = Link('revolute','d', 0, 'a', 0, 'alpha', 0,'modified');
L(2) = Link('revolute','d', 0, 'a', 0.5, 'alpha', 0 ,'modified');
bot = SerialLink(L, 'name', '2D-1-RR');
tool = transl(0.5, 0, 0);
bot.tool = tool;
% bot.teach


q_tool = bot.fkine(q);
q_tool = double(q_tool);
q_tool = squeeze(q_tool(1:3, end, :));

qd_tool = zeros(length(q), 6);
qdd_tool = zeros(length(q), 6);
torques = zeros(length(q), 2);
t_inertia = zeros(length(q), 2);
t_cent = zeros(length(q), 2);
t_cor = zeros(length(q), 2);
t_external = zeros(length(q), 2);
t_grav = zeros(length(q), 2);
J_all = zeros(6, 2, length(q));
val = 9.81;
for i = 1:1:length(q)
    J = jacob(0.5, 0.5, q(i, :));
    J_all(:, :, i) = J;
    qd_tool(i, :) = transpose(J*transpose(qd(i, :)));
    torques(i, :) = subs(Tau, [t1, t1d, t1d2, t1dd, t2, t2d, t2d2, t2dd, g], [q(i, 1), qd(i, 1), qd(i, 1)^2, qdd(i, 1), q(i, 2), qd(i, 2), qd(i, 2)^2, qdd(i, 2), val])';
    t_external(i, :) = subs(Tau, [t1, t1d, t1d2, t1dd, t2, t2d, t2d2, t2dd, g], [q(i, 1), 0, 0, 0, q(i, 2), 0, 0, 0, 0])';
    t_grav(i, :) = subs(Tau, [t1, t1d, t1d2, t1dd, t2, t2d, t2d2, t2dd, g], [q(i, 1), 0, 0, 0, q(i, 2), 0, 0, 0, val])' - t_external(i, :);
    t_inertia(i, :) = subs(Tau, [t1, t1d, t1d2, t1dd, t2, t2d, t2d2, t2dd, g], [q(i, 1), 0, 0, qdd(i, 1), q(i, 2), 0, 0, qdd(i, 2), val])' - t_grav(i, :) - t_external(i, :);
    t_cent(i, :) = subs(Tau, [t1, t1d, t1d2, t1dd, t2, t2d, t2d2, t2dd, g], [q(i, 1), 0, qd(i, 1)^2, 0, q(i, 2), 0, qd(i, 2)^2, 0, val])' - t_grav(i, :) - t_external(i, :);
    t_cor(i, :) = subs(Tau, [t1, t1d, t1d2, t1dd, t2, t2d, t2d2, t2dd, g], [q(i, 1), qd(i, 1), 0, 0, q(i, 2), qd(i, 2), 0, 0, val])' - t_grav(i, :) - t_external(i, :);
end
dt = time(2) - time(1);
dJ = gradient(J_all, dt);

for i = 1:1:length(q)
    qdd_tool(i, :) = transpose(dJ(:, :, i)*transpose(qd(i, :))) + transpose(J_all(:, :, i)*transpose(qdd(i, :)));
end

%%
torques_g = torques;
torques_no_g = torques - t_grav;
torques_no_g./torques_g
figure
hold on
% plot(time, torques)
% plot(time, torques_no_g)
plot(time, torques_no_g./torques_g)
title('Ratio of Joint Torques w/o Gravity to Joint Torques w/ Gravity')
xlabel('Time (s)')
% ylabel('Torque (N*m)')
legend('J1', 'J2')
hold off
%% Plots
close all
% figure
% hold on
% plot(time, q_tool)
% plot(time, qd(:, :))
% title('End Effector Positions Over Time')
% xlabel('Time (s)')
% ylabel('Position (m)')
% legend('X', 'Y', 'Z')
% hold off
% 
% figure
% hold on
% plot(time, qd_tool(:, 1:3))
% title('End Effector Velocity Over Time')
% xlabel('Time (s)')
% ylabel('Velocity (m)')
% legend('X', 'Y', 'Z')
% hold off
% 
% figure
% hold on
% plot(time, qdd_tool(:, 1:3))
% title('End Effector Acceleration Over Time')
% xlabel('Time (s)')
% ylabel('Acceleration (m/s^2)')
% legend('X', 'Y', 'Z')
% hold off

figure
hold on
plot(time, torques)
title('Total Joint Torques Over Time')
xlabel('Time (s)')
ylabel('Torque (N*m)')
legend('J1', 'J2')
hold off

figure
hold on
plot(time, torques(:, 1))
plot(time, t_grav(:, 1))
plot(time, t_inertia(:, 1))
plot(time, t_cor(:, 1))
plot(time, t_cent(:, 1))
plot(time, t_external(:, 1))
title('Joint Torque Components for J1 Over Time')
xlabel('Time (s)')
ylabel('Torque (N*m)')
legend('Total', 'Gravity', 'Inertia', 'Coriolis', 'Centrifugal', 'Applied')
hold off

figure
hold on
plot(time, torques(:, 2))
plot(time, t_grav(:, 2))
plot(time, t_inertia(:, 2))
plot(time, t_cor(:, 2))
plot(time, t_cent(:, 2))
plot(time, t_external(:, 2))
title('Joint Torque Components for J2 Over Time')
xlabel('Time (s)')
ylabel('Torque (N*m)')
legend('Total', 'Gravity', 'Inertia', 'Coriolis', 'Centrifugal', 'Applied')
hold off


%% Plots
figure
hold on
plot(time, q_tool)
% plot(time, qd(:, :))
title('End Effector Positions (No Gravity)')
xlabel('Time (s)')
ylabel('Position (m)')
legend('X', 'Y', 'Z')
hold off

figure
hold on
plot(time, qd_tool(:, 1:3))
title('End Effector Velocity (No Gravity)')
xlabel('Time (s)')
ylabel('Velocity (m)')
legend('X', 'Y', 'Z')
hold off

figure
hold on
plot(time, qdd_tool(:, 1:3))
title('End Effector Acceleration (No Gravity)')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
legend('X', 'Y', 'Z')
hold off

figure
hold on
plot(time, torques)
title('Total Joint Torques (No Gravity)')
xlabel('Time (s)')
ylabel('Torque (N*m)')
legend('J1', 'J2')
hold off

figure
hold on
plot(time, torques(:, 1))
plot(time, t_grav(:, 1))
plot(time, t_inertia(:, 1))
plot(time, t_cor(:, 1))
plot(time, t_cent(:, 1))
title('Joint Torque Components for J1 (No Gravity)')
xlabel('Time (s)')
ylabel('Torque (N*m)')
legend('Total', 'Gravity', 'Inertia', 'Coriolis', 'Centrifugal')
hold off

figure
hold on
plot(time, torques(:, 2))
plot(time, t_grav(:, 2))
plot(time, t_inertia(:, 2))
plot(time, t_cor(:, 2))
plot(time, t_cent(:, 2))
title('Joint Torque Components for J2 (No Gravity)')
xlabel('Time (s)')
ylabel('Torque (N*m)')
legend('Total', 'Gravity', 'Inertia', 'Coriolis', 'Centrifugal')
hold off

%% functions

function J = jacob(L1, L2, joint_angles)
T_0_1 = trotx(0)*transl(0, 0, 0)*trotz(joint_angles(1))*transl(0, 0, 0);
T_1_2 = trotx(0)*transl(L1, 0, 0)*trotz(joint_angles(2))*transl(0, 0, 0);
T_2_tool = trotx(0)*transl(L2, 0, 0)*trotz(0)*transl(0, 0, 0);
T_0_tool = (T_0_1*T_1_2*T_2_tool);

J_exp = zeros(6, 2);
P_0_2 = T_0_tool(1:3, end);
J_exp(1:3, 1) = cross(T_0_1(1:3, end-1), P_0_2 - T_0_1(1:3, end));
T_0_2 = T_0_1*T_1_2;
J_exp(1:3, 2) = cross(T_0_2(1:3, end-1), P_0_2 - T_0_2(1:3, end));
J_exp(4:end, 1) = T_0_1(1:3, end-1);
J_exp(4:end, 2) = T_0_2(1:3, end-1);

% J = J_exp(1:2, :);
J = J_exp;
end

function joint_angles = ik(L1, L2, pt)
% twoR_inverse_kinematics:
%   Computes the inverse kinematics for a 2-link planar arm with link lengths
%   L1 and L2. The end-effector position is (x, y).
%
%   The forward kinematics for joints q1, q2 are:
%       x = L1*cos(q1) + L2*cos(q1 + q2)
%       y = L1*sin(q1) + L2*sin(q1 + q2)
%
%   The solution for q2:
%       q2 = +/- acos( (x^2 + y^2 - L1^2 - L2^2) / (2*L1*L2) )
%
%   Then q1 is found by:
%       q1 = atan2(y, x) - atan2(L2*sin(q2), L1 + L2*cos(q2))
%
%   There can be 0, 1, or 2 feasible solutions, depending on the reachability
%   and geometry. The function returns a 2 x 2 matrix 'sol', where each row is
%   a solution [q1, q2]. If a solution is not feasible, NaN is returned.
%
% Inputs:
%   x, y  - desired end-effector coordinates in the plane
%   L1, L2 - link lengths
%
% Output:
%   sol - a 2x2 matrix of solutions, where:
%             sol(1,:) = [q1_up,   q2_up  ]  (elbow "up" or "down" depending on sign)
%             sol(2,:) = [q1_down, q2_down]
%         If no solution for that branch, the row is [NaN, NaN].

% 1) Compute r^2 = x^2 + y^2
x = pt(1);
y = pt(2);
r2 = x^2 + y^2;

% 2) Check for feasibility
%    The maximum reach is (L1+L2), the minimum reach is |L1-L2|.
%    If r2 > (L1+L2)^2 or r2 < (L1-L2)^2 => no real solutions.
if r2 > (L1 + L2)^2 || r2 < (L1 - L2)^2
    % Entirely unreachable
    joint_angles = [NaN NaN; NaN NaN];
    return;
end

% 3) Compute cos(q2)
c2 = (r2 - L1^2 - L2^2) / (2*L1*L2);

% Numerical issues can make c2 slightly out of [-1, 1] range
if c2 > 1,  c2 = 1;  end
if c2 < -1, c2 = -1; end

% 4) Possible q2 solutions
q2a = acos(c2);   % "elbow" one side
q2b = -acos(c2);  % "elbow" the other side

% 5) Compute corresponding q1 for each q2
%    q1 = atan2(y, x) - atan2(L2*sin(q2), L1 + L2*cos(q2))
%    We'll define a small helper:
q1a = atan2(y, x) - atan2(L2*sin(q2a), L1 + L2*cos(q2a));
q1b = atan2(y, x) - atan2(L2*sin(q2b), L1 + L2*cos(q2b));

% 6) Construct solutions
%    sol(1,:) => solution with q2a
%    sol(2,:) => solution with q2b
% sol = [q1a, q2a;
%        q1b, q2b];
joint_angles = [q1b, q2b];
end