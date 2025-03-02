
clc
close all
clear

L1 = linspace(0, 500, 51);
L2 = linspace(0, 500, 51);

d = 200;

w = linspace(0, 140, 15) + d;
h = linspace(-50, 50, 11);

dw = 15;
dh = 11;

C_max = zeros(1, 3);
for i = 1:1:51
    for j = 1:1:51
        if L1(i)^2 + L2(j)^2 < (d+w(end))^2 + (h(end)/2)^2 || abs(L1(i)-L2(j)) > d
            continue
        else
            % C = optimize(L1(i), L2(j), w, h, dw, dh)
            % if C(1) > C_max(1)
                % C_max = C
            % end
        end
    end
end
optimize(200, 144, w, h, dw, dh)

function C = optimize(L1, L2, w, h, dw, dh)

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
% J_exp(1:3, 3) = cross(T_0_tool(1:3, end-1), P_0_2 - T_0_tool(1:3, end));
J_exp(4:end, 1) = T_0_1(1:3, end-1);
J_exp(4:end, 2) = T_0_2(1:3, end-1);
% J_exp(4:end, 3) = T_0_tool(1:3, end-1);

J_exp = simplify(J_exp);
A = J_exp*transpose(J_exp);


L_1 = Link('revolute', 'd', 0.0, 'a', 0.0, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
L_2 = Link('revolute', 'd', 0.0, 'a', L1, 'alpha', 0.0, 'modified', 'qlim', [-pi, pi]);
L_3 = transl(L2, 0, 0);
bot = SerialLink([L_1, L_2], 'tool', L_3);

k = zeros(1, dw*dh);
for i = 1:1:dh
    for j = 1:1:dw
        joint_angles = ik(L1, L2, [w(j), h(i)]);
        if any(isnan(joint_angles))
            return
        end
        plot = false;
        if plot
            hold on;
            bot.plot(joint_angles, 'workspace', [-200, 400, -75, 400, 0, 0], 'view', [0, 90], 'scale', 0.1)
            patch([200, 340, 340, 200], [50, 50, -50, -50], [0, 0, 0, 0], 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 2);
            plot3(w(j), h(i), 0, 'ro', 'MarkerSize', 4, 'LineWidth', 1);
        end
        lambda = eig(double(subs(A, [t1, t2], joint_angles)));
        % return
        % min_ = min(lambda)
        % max_ = max(lambda)
        % k_ = sqrt(min(lambda)/max(lambda))
        k(j + (i-1)*dh) = sqrt(min(lambda)/max(lambda));
    end
end
k
c = sum(k);
b = min(k);
a = L1^3 + L2^3;

C = [(c*b/a), L1, L2];
end

function joint_angles = ik(L1, L2, pt)
theta = atan2(pt(2), pt(1));
dist = pt(1)^2 + pt(2)^2;
c_t2 = (L1^2 + L2^2 - dist)/(2*L1*L2);
if abs(c_t2) > 1
    joint_angles = NaN;
    return
end
s_t2 = sqrt(1-c_t2^2);
t2 = atan2(s_t2, c_t2) + pi;
% t2_alt = atan2(-s_t2, c_t2);
s_t1 = (L2/sqrt(dist))*s_t2;
c_t1 = sqrt(1-s_t1^2);
t1 = atan2(s_t1, c_t1) + theta;
% t1_alt = atan2(-s_t1, c_t1) + theta;
joint_angles = [t1, t2];
end