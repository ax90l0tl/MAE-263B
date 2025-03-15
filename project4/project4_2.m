
clc
close all
clear
n = 51;
L1 = linspace(0, 500, n);
L2 = linspace(0, 500, n);

d = 200;

w = linspace(0, 140, 15) + d;
h = linspace(-50, 50, 11);

dw = 15;
dh = 11;

C = zeros(n, n);
for i = 1:1:n
    for j = 1:1:n
        C(i, j) = optimize(L1(i), L2(j), w, h, dw, dh);
    end
end
% optimize(.300, .210, w, h, dw, dh);

function C = optimize(L1, L2, w, h, dw, dh)

% syms t1 t2
% T_0_1 = trotx(0)*transl(0, 0, 0)*trotz(t1)*transl(0, 0, 0);
% T_1_2 = trotx(0)*transl(L1, 0, 0)*trotz(t2)*transl(0, 0, 0);
% T_2_tool = trotx(0)*transl(L2, 0, 0)*trotz(0)*transl(0, 0, 0);
% T_0_tool = (T_0_1*T_1_2*T_2_tool);
% 
% J_exp = sym(zeros(6, 2));
% P_0_2 = T_0_tool(1:3, end);
% J_exp(1:3, 1) = cross(T_0_1(1:3, end-1), P_0_2 - T_0_1(1:3, end));
% T_0_2 = T_0_1*T_1_2;
% J_exp(1:3, 2) = cross(T_0_2(1:3, end-1), P_0_2 - T_0_2(1:3, end));
% J_exp(4:end, 1) = T_0_1(1:3, end-1);
% J_exp(4:end, 2) = T_0_2(1:3, end-1);
% 
% J_exp = simplify(J_exp);
% J_square = J_exp(1:2, :);
% simplify(det(J_square))
% A = J_square*transpose(J_square);
% simplify(A)

% L_1 = Link('revolute', 'd', 0.0, 'a', 0.0, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
% L_2 = Link('revolute', 'd', 0.0, 'a', L1, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
% L_3 = transl(L2, 0, 0);
% bot = SerialLink([L_1, L_2], 'tool', L_3);

k = zeros(dh, dw);
for i = 1:1:dh
    for j = 1:1:dw
        joint_angles = ik(L1, L2, [w(j), h(i)]);
        J_square = jacob(L1, L2, joint_angles);
        A = J_square*transpose(J_square);
        if any(isnan(joint_angles))
            C = 0;
            return
        end

        % hold on;
        % bot.plot(joint_angles, 'workspace', [-200, 400, -75, 400, 0, 0], 'view', [0, 90], 'scale', 0.1)
        % patch([200, 340, 340, 200], [50, 50, -50, -50], [0, 0, 0, 0], 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 2);
        % plot3(w(j), h(i), 0, 'ro', 'MarkerSize', 4, 'LineWidth', 1);

        lambda = eig(A);
        k(i, j) = sqrt(min(abs(lambda))/max(abs(lambda)));
    end
end
c = sum(k, 'all');
b = min(k, [], 'all');
a = L1^3 + L2^3;

C = (c*b/a);
end

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

J = J_exp(1:2, :);
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

%%
surf(L1, L2, C);
xlabel('L1 (mm)')
ylabel('L2 (mm)')
title('Goal Value')
colorbar

%%
close all
clc
m = max(C, [], 'all');

[j, k] = find(C == m);

L1_best = L1(j);
L2_best = L2(k);

nonzero_vals = C(C > 0);
[min_val, idx] = min(nonzero_vals);
[j, k] = find(C == min_val, 1);

L1_worst = L1(j);
L2_worst = L2(k);

[goal_best, k_best] = optimize2(L1_best, L2_best, w, h, dw, dh);
[goal_worst, k_worst] = optimize2(L1_worst, L2_worst, w, h, dw, dh);
figure
surf(w, h, k_best);
xlabel('X (mm)')
ylabel('Y (mm)')
title('Best Configuration')
colorbar

figure
surf(w, h, k_worst);
xlabel('X (mm)')
ylabel('Y (mm)')
title('Worst Configuration')
colorbar

function [C, k] = optimize2(L1, L2, w, h, dw, dh)

syms t1 t2
T_0_1 = trotx(0)*transl(0, 0, 0)*trotz(t1)*transl(0, 0, 0);
T_1_2 = trotx(0)*transl(L1, 0, 0)*trotz(t2)*transl(0, 0, 0);
T_2_tool = trotx(0)*transl(L2, 0, 0)*trotz(0)*transl(0, 0, 0);
T_0_tool = (T_0_1*T_1_2*T_2_tool);

J_exp = sym(zeros(6, 2));
P_0_2 = T_0_tool(1:3, end);
J_exp(1:3, 1) = cross(T_0_1(1:3, end-1), P_0_2 - T_0_1(1:3, end));
T_0_2 = T_0_1*T_1_2;
J_exp(1:3, 2) = cross(T_0_2(1:3, end-1), P_0_2 - T_0_2(1:3, end));
J_exp(4:end, 1) = T_0_1(1:3, end-1);
J_exp(4:end, 2) = T_0_2(1:3, end-1);

J_exp = simplify(J_exp);
J_square = J_exp(1:2, :);
% simplify(det(J_square))
A = J_square*transpose(J_square);
% simplify(A)

% L_1 = Link('revolute', 'd', 0.0, 'a', 0.0, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
% L_2 = Link('revolute', 'd', 0.0, 'a', L1, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
% L_3 = transl(L2, 0, 0);
% bot = SerialLink([L_1, L_2], 'tool', L_3);

k = zeros(dh, dw);
min_lamba = 100000000;
for i = 1:1:dh
    for j = 1:1:dw
        joint_angles = ik(L1, L2, [w(j), h(i)]);
        if any(isnan(joint_angles))
            C = 0;
            return
        end

        % hold on;
        % bot.plot(joint_angles, 'workspace', [-200, 400, -75, 400, 0, 0], 'view', [0, 90], 'scale', 0.1)
        % patch([200, 340, 340, 200], [50, 50, -50, -50], [0, 0, 0, 0], 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 2);
        % plot3(w(j), h(i), 0, 'ro', 'MarkerSize', 4, 'LineWidth', 1);

        lambda = eig(double(subs(A, [t1, t2], joint_angles)));
        if min(lambda) < min_lamba
            min_lamba = min(lambda)
            w(j), h(i)
        end
        k(i, j) = sqrt(min((lambda))/max((lambda)));
    end
end
c = sum(k, 'all');
b = min(k, [], 'all');
a = L1^3 + L2^3;

C = (c*b/a);
end