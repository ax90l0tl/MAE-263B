%% A
clc

rho = 2710;
m_cube = rho*0.1^3;
I_cube = (1/12)*m_cube*diag([0.1^2 + 0.1^2, 0.1^2 + 0.1^2, 0.1^2 + 0.1^2])
m_cylinder = rho*((0.02/2)^2)*pi*0.1;
I_cylinder = m_cylinder*diag([(1/12)*(3*(0.02/2)^2 + 0.1^2), (1/12)*(3*(0.02/2)^2 + 0.1^2), 0.5*(0.02/2)^2])
m = m_cube - m_cylinder;
I_A = I_cube - I_cylinder

P_A_CM = [-0.4;0;0];

I_A_CM = I_A + m*(transpose(P_A_CM)*P_A_CM*eye(3) - P_A_CM'.*P_A_CM)

%% B
clc
r = 0.05/2;
l = 0.8 - 0.1;
m_cylinder = rho*(r^2)*pi*0.1;
I_cylinder = m_cylinder*diag([0.5*(r)^2, (1/12)*(3*(r)^2 + l^2), (1/12)*(3*(r)^2 + l^2)]);
I_B_CM = I_cylinder
%% C
clc
m_cube = rho*(0.1^3);
r = 0.05/2;
l = 0.1;
I_cube = (1/12)*m_cube*diag([0.1^2 + 0.1^2, 0.1^2 + 0.1^2, 0.1^2 + 0.1^2]);
m_cylinder = rho*(r^2)*pi*0.1;
I_cylinder = m_cylinder*diag([(1/12)*(3*(r)^2 + l^2), (1/12)*(3*(r)^2 + l^2), 0.5*(r)^2]);
m = m_cube - m_cylinder;
I_C = I_cube - I_cylinder;

P_C_CM = [0.4;0;0];
rot = trotx(pi/4);
rot = rot(1:3, 1:3);
syms t
trotx(t)
I_C_CM = rot*I_C*transpose(rot) + m*(transpose(P_C_CM)*P_C_CM*eye(3) - P_C_CM'.*P_C_CM)

%%
clc
I_total = (I_A_CM + I_B_CM + I_C_CM)