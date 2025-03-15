clc
close all
clear

% ellipsoid(0.25, 0.75)
% ellipsoid(0.5, 0.5);
% ellipsoid(0.75, 0.25)

function C = ellipsoid(L1, L2)

syms t1 t2
T_0_1 = trotx(0)*transl(0, 0, 0)*trotz(t1)*transl(0, 0, 0);
T_1_2 = trotx(0)*transl(L1, 0, 0)*trotz(t2)*transl(0, 0, 0);
T_2_tool = trotx(0)*transl(L2, 0, 0)*trotz(0)*transl(0, 0, 0);
T_0_tool = simplify(T_0_1*T_1_2*T_2_tool);

J_exp = sym(zeros(6, 2));
P_0_2 = T_0_tool(1:3, end);
J_exp(1:3, 1) = cross(T_0_1(1:3, end-1), P_0_2 - T_0_1(1:3, end));
T_0_2 = T_0_1*T_1_2;
J_exp(1:3, 2) = cross(T_0_2(1:3, end-1), P_0_2 - T_0_2(1:3, end));
J_exp(4:end, 1) = T_0_1(1:3, end-1);
J_exp(4:end, 2) = T_0_2(1:3, end-1);

J_exp = simplify(J_exp);
J_square = J_exp(1:2, :);
A = J_square*transpose(J_square);
simplify(A);

% L_1 = Link('revolute', 'd', 0.0, 'a', 0.0, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
% L_2 = Link('revolute', 'd', 0.0, 'a', L1, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
% L_3 = transl(L2, 0, 0);
% bot = SerialLink([L_1, L_2], 'tool', L_3);

x = linspace(0, 1.0, 11);
vel_legend = string([]);
force_legend = string([]);
for i = 1:1:11
    joint_angles = ik(L1, L2, [x(i), 0]);
    if any(isnan(joint_angles))
        continue
    end
    figure(1);
    hold on;
    title("Velocity Ellipsoid")
    % xlim([-3.0, 3.0]);
    % ylim([-1.0, 1.0]);
    % bot.plot(joint_angles, 'view', [0, 90], 'scale', 0.1)
    % bot.fkine(joint_angles);
    plot3(x(i), 0, 0, 'ro', 'MarkerSize', 4, 'LineWidth', 1);
    xlabel('X (m)')
    % velocity ellipse
    [vector, lambda] = eig(subs(A, [t1, t2], joint_angles));
    vector = vector./vecnorm(double(vector), 2, 1);
    lambda = sqrt(diag(lambda));
    manip = [max(lambda)/min(lambda), x(i)]
    theta = linspace(0, 2*pi, 100);
    pos = lambda(1)*vector(:, 1)*cos(theta)+lambda(2)*vector(:, 2)*sin(theta) + [x(i); 0];
    plot3(pos(1, :), pos(2, :), zeros(100, 1), 'LineWidth', 1);
    vel_legend = [vel_legend, "", string(x(i))];
    legend(vel_legend)
    hold off

    figure(2)
    hold on
    title("Force Ellipsoid")
    % xlim([-3.0, 3.0]);
    % ylim([-1.0, 1.0]);
    % bot.plot(joint_angles, 'view', [0, 90], 'scale', 0.1)
    % bot.fkine(joint_angles);
    xlabel('X (m)')
    plot3(x(i), 0, 0, 'ro', 'MarkerSize', 4, 'LineWidth', 1);
    force_legend = [force_legend, ""];
    legend(force_legend)
    % force ellipse
    if joint_angles(2) ~= 0 && joint_angles(2)/pi ~= 1
        [vector, lambda] = eig(inv(subs(transpose(A), [t1, t2], joint_angles)));
        vector = vector./vecnorm(double(vector), 2, 1);
        lambda = sqrt(diag(lambda));
        theta = linspace(0, 2*pi, 100);
        pos = lambda(1)*vector(:, 1)*cos(theta)+lambda(2)*vector(:, 2)*sin(theta) + [x(i); 0];
        plot3(pos(1, :), pos(2, :), zeros(100, 1), 'LineWidth', 1);
        force_legend = [force_legend, string(x(i))];
        legend(force_legend)
        hold off
    end
end
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
    joint_angles = [NaN NaN];
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
joint_angles = [q1a, q2a];
end

