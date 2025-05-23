%% explicit method
clc
clear
close all

syms d1 d2 d3 d4 d5 a1 a2 a3 a4 alpha1 alpha2 alpha3 alpha4
syms t1 t2 t3 d4
% test
% t1 = 0;
% t2 = 0;
% t3 = 0;
% d4 = 0;
% T_0_1 = trotx(0)*transl(0, 0, 0)*trotz(t1)*transl(0, 0, 0.400);
% T_1_2 = trotx(0)*transl(0.325, 0, 0)*trotz(t2)*transl(0, 0, 0);
% T_2_3 = trotx(0)*transl(0.225, 0, 0)*trotz(t3)*transl(0, 0, 0);
% T_3_4 = trotx(0)*transl(0, 0, 0)*trotz(0)*transl(0, 0, d4);
% T_4_tool = trotx(0)*transl(0, 0, 0)*trotz(0)*transl(0, 0, -0.03);
% T_0_tool = (T_0_1*T_1_2*T_2_3*T_3_4*T_4_tool);

% frames
T_0_1 = trotx(0)*transl(0, 0, 0)*trotz(t1)*transl(0, 0, d1);
T_1_2 = trotx(0)*transl(a1, 0, 0)*trotz(t2)*transl(0, 0, 0);
T_2_3 = trotx(0)*transl(a2, 0, 0)*trotz(t3)*transl(0, 0, 0);
T_3_4 = trotx(0)*transl(0, 0, 0)*trotz(0)*transl(0, 0, d4);
T_4_tool = trotx(0)* transl(a4, 0, 0)*trotz(0)*transl(0, 0, d5);
T_0_4 = simplify(T_0_1*T_1_2*T_2_3*T_3_4);

T_all(:,:,1) = T_0_1;
T_all(:,:,2) = T_1_2;
T_all(:,:,3) = T_2_3;
T_all(:,:,4) = T_3_4;
T_all(:,:,5) = T_0_4;
size(T_all);

J_exp = sym(zeros(6, 4));
P_0_4 = T_0_4(1:3, end);
J_exp(1:3, 1) = cross(T_all(1:3, end-1, 1), P_0_4 - T_0_1(1:3, end));
T_0_2 = T_0_1*T_1_2;
T_0_3 = T_0_1*T_1_2*T_2_3;
J_exp(1:3, 2) = cross(T_all(1:3, end-1, 2), P_0_4 - T_0_2(1:3, end));
J_exp(1:3, 3) = cross(T_all(1:3, end-1, 3), P_0_4 - T_0_3(1:3, end));
J_exp(1:3, 4) = T_all(1:3, end-1, 4);
J_exp(4:end, 1) = T_0_1(1:3, 3);
J_exp(4:end, 2) = T_0_2(1:3, 3);
J_exp(4:end, 3) = T_0_3(1:3, 3);

J_exp = simplify(J_exp);
pretty(J_exp)
J_square = J_exp([1:3, end], :);
sing = simplify(det(J_square))
%% Velocity propagation
syms t1_dot t2_dot t3_dot d4_dot
clc
% angular velocity
w_0 = [0;0;0;];
w_1 = T_all(1:3, 1:3, 1)*w_0 + [0;0;t1_dot];
w_2 = T_all(1:3, 1:3, 2)*w_1 + [0;0;t2_dot];
w_3 = T_all(1:3, 1:3, 3)*w_2 + [0;0;t3_dot];
w_4 = T_all(1:3, 1:3, 4)*w_3 + [0;0;0];
% linear velocity
v_0 = [0;0;0;];
v_1 = transpose(T_all(1:3, 1:3, 1))*(cross(w_0, T_all(1:3, end, 1)) + v_0);
v_2 = transpose(T_all(1:3, 1:3, 2))*(cross(w_1, T_all(1:3, end, 2)) + v_1);
v_3 = transpose(T_all(1:3, 1:3, 3))*(cross(w_2, T_all(1:3, end, 3)) + v_2);
v_4 = transpose(T_all(1:3, 1:3, 4))*(cross(w_3, T_all(1:3, end, 4)) + v_3) + [0;0;d4_dot];
% transform to base frame
w_4 = T_all(1:3, 1:3, 5)*w_4;
v_4 = simplify(T_all(1:3, 1:3, 5)*v_4);

all_v = [v_4; w_4];
q_dot = [t1_dot; t2_dot; t3_dot; d4_dot];

J_vel = simplify(jacobian(all_v, q_dot));
pretty(J_vel)
size(J_vel);

J_square = J_vel([1:3, end], :);
sing = simplify(det(J_square))
%% Force Propagation
syms fx fy fz tx ty tz torque1 torque2 torque3 force4
clc

% End effector forces
f_4 = [fx;fy;fz];
t_4 = [tx;ty;tz];
% Forces
f_3 = simplify((T_3_4(1:3, 1:3))*f_4);
f_2 = simplify((T_2_3(1:3, 1:3))*f_3);
f_1 = simplify((T_1_2(1:3, 1:3))*f_2);
% Torques
t_3 = simplify((T_3_4(1:3, 1:3))*t_4 + cross(T_3_4(1:3, end), f_3));
t_2 = simplify((T_2_3(1:3, 1:3))*t_3 + cross(T_2_3(1:3, end), f_2));
t_1 = simplify((T_1_2(1:3, 1:3))*t_2 + cross(T_1_2(1:3, end), f_1));
% joint torques
jt_1 = simplify(collect(transpose(t_1)*[0;0;1], [fx, fy, fz, tx, ty, tz]));
jt_2 = simplify(collect(transpose(t_2)*[0;0;1], [fx, fy, fz, tx, ty, tz]));
jt_3 = simplify(collect(transpose(t_3)*[0;0;1], [fx, fy, fz, tx, ty, tz]));
jt_4 = simplify(collect(transpose(f_4)*[0;0;1], [fx, fy, fz, tx, ty, tz]));

wrench = [fx;fy;fz;tx;ty;tz];
% building jacobian
tq_1 = equationsToMatrix(jt_1, wrench);
tq_2 = equationsToMatrix(jt_2, wrench);
tq_3 = equationsToMatrix(jt_3, wrench);
tq_4 = equationsToMatrix(jt_4, wrench);

J_force_T = [tq_1;tq_2;tq_3;tq_4];
% Jacobian in end effector frame
J_force = transpose(J_force_T);


transform = [(T_0_4(1:3, 1:3)), zeros(3,3);
    zeros(3,3),(T_0_4(1:3,1:3))];
% Jacobian in base frame
J_force_0 = simplify(transform*J_force)

J_square = J_force_0([1:3, end], :)
sing = simplify(det(J_square))