clc; clear all; close all;
sympref('FloatingPointOutput',true);
format short
%% Construct Arm
syms q1(t) q2(t) % time varying joint parameters
q = [q1(t) , q2(t)];
syms q1d(t) q2d(t) 
qd = [q1d(t) , q2d(t)];
syms q1dd(t) q2dd(t)
qdd = [q1dd(t) , q2dd(t)];


thickness = 0.005; density = 2700; % m | kg/m^3
L_val = [0.5 , 0.5]; % m
syms L1 L2
L = [L1 , L2];
OD = [0.1 , 0.1]; % m 
ID = OD - thickness; % m 
I_x = sym(zeros(1,length(L))); I_y = sym(zeros(1,length(L))); I_z = sym(zeros(1,length(L))); 
V = sym(zeros(1,length(L))); M = sym(zeros(1,length(L)));
syms g
% g = [ 0 ; -9.81 ; 0 ; ];

for i=1:length(L)
    V(i) = pi * (OD(i)^2 - ID(i)^2) * L(i);
    M(i) = V(i) * density;
    I_x(i) = 0.5 * density * V(i) * ( OD(i)^2 + ID(i)^2 ); 
    I_y(i) = (1/12) * density * V(i) * ( 3 * ( OD(i)^2 + ID(i)^2 ) + L(i)^2 );
    I_z(i) = I_y(i);
end

% DH Parameters and FK
DH_table = [ 0 , q1(t) , 0    , 0         ;
             0 , q2(t) , L(1) , 0         ;
             0 , 0         , L(2) , 0         ;];

T.T01 = DH_Matrix(DH_table(1,:));
T.T12 = DH_Matrix(DH_table(2,:));
T.T23 = DH_Matrix(DH_table(3,:));

T.T02 = simplify(T.T01 * T.T12);
T.T03 = simplify(T.T01 * T.T12 * T.T23);
T_end = T.T03; 
P1_function = @(Q) subs(T.T01(1:3,4), q, Q');
P2_function = @(Q) subs(T.T02(1:3,4), q, Q');
P3_function = @(Q) subs(T.T03(1:3,4), q, Q');

% Hard-Coding Stuff for Sanity
R01 = T.T01(1:3,1:3);
R12 = T.T12(1:3,1:3);
R23 = T.T23(1:3,1:3);
R10 = R01.';
R21 = R12.';
R32 = R23.';
R03 = T.T03(1:3,1:3);
R30 = R03.';
P01 = T.T01(1:3,4);
P12 = T.T12(1:3,4);
P23 = T.T23(1:3,4);
%% Newton-Euler
syms g
syms m1 m2
m = [m1 , m2];
% syms L1 L2
% L = [L1 , L2];
IC1 = (1/12) * m1 * L1^2 * [ 0 0 0 ; 0 1 0 ; 0 0 1 ];
IC2 = (1/12) * m2 * L2^2 * [ 0 0 0 ; 0 1 0 ; 0 0 1 ];
syms Lc1 Lc2
Lcom = [Lc1 , Lc2];


omega = sym(zeros(3,length(L)+1)); omegad = sym(zeros(3,length(L)+1));
vel_d = sym(zeros(3,length(L)+1));   vel_d_com = sym(zeros(3,length(L)+1));
omega(:,1) = [ 0; 0; 0;];     omegad(:,1) = [ 0; 0; 0;]; 
vel_d(:,1) = [0;g;0];                vel_d_com(:,1) = [ 0; 0; 0;]; 
F = sym(zeros(3,length(L)));  N = sym(zeros(3,length(L))); 

%Forward Iterations
for i=1:length(L)
    index = strcat('T',num2str(i-1),num2str(i));
    next_index = strcat('T',num2str(i),num2str(i+1));
    R = T.(index)(1:3,1:3);
    P = T.(index)(1:3,4);
    P_com = [Lcom(i);0;0];
    omega(:,i+1) = simplify(R.' * omega(:,i) + qd(i)*R(:,3));
    omegad(:,i+1) = simplify(R.' * omegad(:,i) + cross( R.' * omega(:,i) , qd(i) * R(:,3) )  + qdd(i) * R(:,3));
    vel_d(:,i+1) = simplify(R.' * (cross( omegad(:,i),P) + cross( omega(:,i) , cross(omega(:,i),P)) + vel_d(:,i)));
    vel_d_com(:,i+1) = cross(omegad(:,i+1), P_com) + cross( omega(:,i+1) , cross( omega(:,i+1), P_com)) + vel_d(:,i+1); 
    F(:,i) = m(i) * vel_d_com(:,i+1);
    I = (1/12) * m(i) * L(i)^2 * [ 0 0 0 ; 0 1 0 ; 0 0 1 ];
    N(:,i) =  I * omegad(:,i+1) + cross(omegad(:,i+1),I*omegad(:,i+1));
end

f = sym(zeros(3,length(L)+1)); n = sym(zeros(3,length(L)+1));
f_app = [-10;0;0]; n_app = [0;0;10];
f(:,end) = R30*f_app; n(:,end) = R30*n_app; % Zero Initial Forces for

tau = sym(zeros(length(L),1));
% Backwards
for i=2:-1:1
    index = strcat('T',num2str(i),num2str(i+1));
    R = T.(index)(1:3,1:3);
    P = T.(index)(1:3,4);
    P_com = [Lcom(i);0;0]; %0.5 * T.(index)(1:3,4);
    f(:,i) = R * f(:,i+1) + F(:,i);
    n(:,i) = N(:,i) + R*n(:,i+1) + cross(P_com,F(:,i)) +  cross(P,R*f(:,i+1));
    tau(i) = simplify(transpose(n(:,i))*[0;0;1]);
end

%% Lagrange Formulation
syms m1 m2
m = [m1 , m2];
syms Iz1 Iz2
Iz = [Iz1, Iz2];
syms g
syms Lc1 Lc2

T0COM1 = T.T01*transl([Lc1,0,0]);
T0COM2 = T.T02*transl([Lc2,0,0]);
P0COM1 = T0COM1(1:3,4);
P0COM2 = T0COM2(1:3,4);


V0COM1 = subs(diff(P0COM1,t),diff(q,t),qd);
W0COM1 = q1d(t);
K1 = simplify(0.5 * m(1) * (V0COM1.'*V0COM1) + 0.5 * Iz(1) * W0COM1^2);
P1 = simplify(m(1) * g * P0COM1(2));


V0COM2 = subs(diff(P0COM2,t),diff(q,t),qd);
W0COM2 = q1d(t)+q2d(t);

K2 = simplify(0.5 * m(2) * (V0COM2.'*V0COM2) + 0.5 * Iz(2) * W0COM2^2);
P2 = simplify(m(2) * g * P0COM2(2));

K = K1+K2;
P = (P1 + P2);
LG = K - P;
EOM = [ subs(diff(diff(K,qd(1)),t),[diff(qd,t),diff(q,t)],[qdd,qd]) - diff(K,q(1)) + diff(P,q(1)); 
        subs(diff(diff(K,qd(2)),t),[diff(qd,t),diff(q,t)],[qdd,qd]) - diff(K,q(2)) + diff(P,q(2));    ];

J_ext = [ -L(1)*sin(q(1)) - L(2)*sin(q(1)+q(2)) , -L(2)*sin(q(1)+q(2));
           L(1)*cos(q(1)) + L(2)*cos(q(1)+q(2)) ,  L(2)*cos(q(1)+q(2));
           1                                    ,  1                  ;];

F_applied = f_app+n_app; % A force applied in -X and Positive Moment Applied
EOM_ext = EOM + J_ext.'*(F_applied);
%% Substitute Back in Values for Both EOM
M = subs(M,L,L_val);
Lcom_val = 0.5*L_val;
Iz_val = subs([IC1(3,3),IC2(3,3)],[m,L],[M,L_val]);
EOM_val = subs(EOM_ext,[m,g,L,Lcom,Iz],[M,9.81,L_val,Lcom_val,Iz_val]); 

tau_val = subs(tau,[g,L,Lcom,m],[9.81,L_val,Lcom_val,M]);
% Comparison
if simplify(tau_val-EOM_val)~=0
    fprintf('\n Equations of Motion do not Match!')
else
    fprintf('\n Equations of Motion Match!')
end

% Hard Code the Expressions
inertial = [ 0.0833*L1^2*m1*q1dd(t) + L1^2*m2*q1dd(t) + 0.0833*L2^2*m2*q1dd(t) + 0.0833*L2^2*m2*q2dd(t) + Lc1^2*m1*q1dd(t) + Lc2^2*m2*q1dd(t) + Lc2^2*m2*q2dd(t) +  2*L1*Lc2*m2*cos(q2(t))*q1dd(t) + L1*Lc2*m2*cos(q2(t))*q2dd(t);
    0.0833*L2^2*m2*q1dd(t) + 0.0833*L2^2*m2*q2dd(t) + Lc2^2*m2*q1dd(t) + Lc2^2*m2*q2dd(t) + L1*Lc2*m2*cos(q2(t))*q1dd(t) ];
coriolis = [ - 2*L1*Lc2*m2*sin(q2(t))*q1d(t)*q2d(t) ; 0 ];
centrifugal = [ - L1*Lc2*m2*sin(q2(t))*q2d(t)^2 ; 
    L1*Lc2*m2*sin(q2(t))*q1d(t)^2];
gravitational = [ L1*g*m2*cos(q1(t)) + Lc1*g*m1*cos(q1(t)) + Lc2*g*m2*cos(q1(t))*cos(q2(t)) - Lc2*g*m2*sin(q1(t))*sin(q2(t)); 
    Lc2*g*m2*cos(q1(t))*cos(q2(t)) - Lc2*g*m2*sin(q1(t))*sin(q2(t)) ];

inertial_val_g = subs(inertial,[g,L,Lcom,m],[9.81,L_val,Lcom_val,M]);
coriolis_val_g = subs(coriolis,[g,L,Lcom,m],[9.81,L_val,Lcom_val,M]);
centrifugal_val_g = subs(centrifugal,[g,L,Lcom,m],[9.81,L_val,Lcom_val,M]);
gravitational_val_g = subs(gravitational,[g,L,Lcom,m],[9.81,L_val,Lcom_val,M]);

inertial_val_no_g = subs(inertial,[g,L,Lcom,m],[0,L_val,Lcom_val,M]);
coriolis_val_no_g = subs(coriolis,[g,L,Lcom,m],[0,L_val,Lcom_val,M]);
centrifugal_val_no_g = subs(centrifugal,[g,L,Lcom,m],[0,L_val,Lcom_val,M]);
gravitational_val_no_g = subs(gravitational,[g,L,Lcom,m],[0,L_val,Lcom_val,M]);
1

%% Joint Space Trajectory Generation - Gravity
time = 0:0.1:4;
Q_start = fastIK([0.1,0],L_val);
Q_end = fastIK([0.9,0],L_val);
[Q,QD,QDD,QPP] = quinticpolytraj([Q_start,Q_end], [0 , 4], time );

tau_applied = zeros(2,length(time));
tau_inertial = zeros(2,length(time));
tau_coriolis = zeros(2,length(time));
tau_centrifugal =  zeros(2,length(time));
tau_gravitational = zeros(2,length(time));
tau_total = zeros(2,length(time));

tau_inertial_no_g = zeros(2,length(time));
tau_coriolis_no_g = zeros(2,length(time));
tau_centrifugal_no_g =  zeros(2,length(time));
tau_gravitational_no_g = zeros(2,length(time));
tau_total_no_g = zeros(2,length(time));

position = zeros(2,length(time));
velocity = zeros(2,length(time));
acceleration = zeros(2,length(time));

for i=1:length(time)
    % End Effector Task Space Kinematics
    position(:,i) = fastP(Q(:,i).',L_val);
    velocity(:,i) = fastJv(Q(:,i).',L_val) * QD(:,i);
    acceleration(:,i) = fastJv(Q(:,i).',L_val) * QDD(:,i) + fastJv_dot(Q(:,i).',QD(:,i).',L_val) * QD(:,i);

    % Torque Calculation 
    tau_applied(:,i) = fastJ(Q(:,i),L_val).'*F_applied;
    tau_inertial(:,i) = subs(inertial_val_g,[q qd qdd],[Q(:,i).',QD(:,i).',QDD(:,i).'] );
    tau_coriolis(:,i) = subs(coriolis_val_g,[q qd qdd],[Q(:,i).',QD(:,i).',QDD(:,i).'] );
    tau_centrifugal(:,i) = subs(centrifugal_val_g,[q qd qdd],[Q(:,i).',QD(:,i).',QDD(:,i).'] );
    tau_gravitational(:,i) = subs(gravitational_val_g,[q qd qdd],[Q(:,i).',QD(:,i).',QDD(:,i).'] );
    tau_total(:,i) = tau_applied(:,i) + tau_inertial(:,i) + tau_coriolis(:,i) + tau_centrifugal(:,i) + tau_gravitational(:,i);

    tau_inertial_no_g(:,i) = subs(inertial_val_no_g,[q qd qdd],[Q(:,i).',QD(:,i).',QDD(:,i).'] );
    tau_coriolis_no_g(:,i) = subs(coriolis_val_no_g,[q qd qdd],[Q(:,i).',QD(:,i).',QDD(:,i).'] );
    tau_centrifugal_no_g(:,i) = subs(centrifugal_val_no_g,[q qd qdd],[Q(:,i).',QD(:,i).',QDD(:,i).'] );
    tau_gravitational_no_g(:,i) = subs(gravitational_val_no_g,[q qd qdd],[Q(:,i).',QD(:,i).',QDD(:,i).'] );
    tau_total_no_g(:,i) = tau_applied(:,i) + tau_inertial_no_g(:,i) + tau_coriolis_no_g(:,i) + tau_centrifugal_no_g(:,i) + tau_gravitational_no_g(:,i);
end
%% Plotting
close all

figure(1)
hold on
plot(time,position(1,:))
plot(time,velocity(1,:))
plot(time,acceleration(1,:))
legend('Position [m]','Velocity [m/s]','Acceleration [m/s^2]')
xlabel('Time [s]')
ylabel('Value')
title('Task Space Kinematics of End Effector in X-Axis')

figure(2)
hold on
plot(time,tau_total(1,:))
plot(time,tau_inertial(1,:))
plot(time,tau_coriolis(1,:))
plot(time,tau_centrifugal(1,:))
plot(time,tau_gravitational(1,:))
plot(time,tau_applied(1,:))
lgd = legend('\tau_{total}','\tau_{inertial}','\tau_{coriolis}', '\tau_{centrifugal}','\tau_{gravitational}','\tau_{applied}');
lgd.FontSize = 14;
xlabel('Time [s]')
ylabel('Torque [N-m]')
title('Joint Torque About Motor 1 with Gravity')

figure(3)
hold on
plot(time,tau_total(2,:))
plot(time,tau_inertial(2,:))
plot(time,tau_coriolis(2,:))
plot(time,tau_centrifugal(2,:))
plot(time,tau_gravitational(2,:))
plot(time,tau_applied(2,:))
lgd = legend('\tau_{total}','\tau_{inertial}','\tau_{coriolis}', '\tau_{centrifugal}','\tau_{gravitational}','\tau_{applied}');
lgd.FontSize = 14;
xlabel('Time [s]')
ylabel('Torque [N-m]')
title('Joint Torque About Motor 2 with Gravity')

figure(4)
hold on 
plot(time,tau_total_no_g(1,:))
plot(time,tau_inertial_no_g(1,:))
plot(time,tau_coriolis_no_g(1,:))
plot(time,tau_centrifugal_no_g(1,:))
plot(time,tau_gravitational_no_g(1,:))
plot(time,tau_applied(1,:))
lgd = legend('\tau_{total}','\tau_{inertial}','\tau_{coriolis}', '\tau_{centrifugal}','\tau_{gravitational}','\tau_{applied}');
lgd.FontSize = 14;
xlabel('Time [s]')
ylabel('Torque [N-m]')
title('Joint Torque About Motor 1 without Gravity')

figure(5)
hold on
plot(time,tau_total_no_g(2,:))
plot(time,tau_inertial_no_g(2,:))
plot(time,tau_coriolis_no_g(2,:))
plot(time,tau_centrifugal_no_g(2,:))
plot(time,tau_gravitational_no_g(2,:))
plot(time,tau_applied(2,:))
lgd = legend('\tau_{total}','\tau_{inertial}','\tau_{coriolis}', '\tau_{centrifugal}','\tau_{gravitational}','\tau_{applied}');
lgd.FontSize = 14;
xlabel('Time [s]')
ylabel('Torque [N-m]')
title('Joint Torque About Motor 2 without Gravity')

figure(6)
hold on
plot(time,tau_total_no_g(1,:)./tau_total(1,:))
plot(time,tau_total_no_g(2,:)./tau_total(2,:))
legend('Ratio of \tau Without Gravity to Gravity for Motor 1','Ratio of \tau Without Gravity to Gravity for Motor 2')
xlabel('Time [s]')
ylabel('Ratio')
title('Comparisons of Total Joint Torques With and Without Gravity')

%% Problem 2
rho = 2700; %kg/m^3
L_cyl = 0.80;
D_cyl = 0.05;
r_cyl = D_cyl/2;
L_prism = 0.1;
W_prism = 0.1;
H_prism = 0.1;
D_hole = 0.02;
r_hole = D_hole/2;
L_hole = 0.1;

% Cylinder
V_cyl = L_cyl*pi*r_cyl^2;
Ix_cyl = 0.5 * V_cyl * rho * r_cyl^2;
Iy_cyl = (1/12) * V_cyl * rho * (3*r_cyl^2 + L_cyl^2 );
Iz_cyl = Iy_cyl;
I_cyl = diag([Ix_cyl,Iy_cyl,Iz_cyl]);

% Prism
V_prism = L_prism * W_prism * H_prism;
Ix_prism = 1/12 * V_prism * rho * (L_prism^2 + H_prism^2);
Iy_prism = 1/12 * V_prism * rho * (W_prism^2 + H_prism^2);
Iz_prism = 1/12 * V_prism * rho * (L_prism^2 + W_prism^2);
I_prism = diag([Ix_prism, Iy_prism, Iz_prism]);

% Hole 
V_hole = L_hole * pi * r_hole^2;
Ix_hole = (1/12) * V_hole * rho * (3*r_hole^2 + L_hole^2 );
Iy_hole = Ix_hole;
Iz_hole = 0.5 * V_hole * rho * r_hole^2;
I_hole = diag([Ix_hole, Iy_hole, Iz_hole]);

I_box = I_prism - I_hole;

r_a = [-0.4 ; 0 ; 0 ];
I_box_a = I_box + rho*(V_prism - V_hole) * (r_a.'*r_a*eye(3,3) - r_a*r_a.');

r_c = [ 0.4 ; 0 ; 0 ];
Rc = rotx(-45,"deg");
I_box_c_rot = Rc*I_box*Rc.';
I_box_c = I_box_c_rot + rho*(V_prism - V_hole) * (r_c.'*r_c*eye(3,3) - r_c*r_c.');

I_tot = I_box_a + I_cyl + I_box_c

%% Defunct Testing Code
% PC1 = [L(1)/2;0;0];
% PC2 = [L(2)/2;0;0];
% 
% IC1 = (1/12) * m(1) * L(1)^2 * [ 0 0 0 ; 0 1 0 ; 0 0 1 ];
% IC2 = (1/12) * m(2) * L(2)^2 * [ 0 0 0 ; 0 1 0 ; 0 0 1 ];
% 
% f3 = zeros(3,1);
% n3 = zeros(3,1);
% 
% w0 = zeros(3,1);
% wd0 = zeros(3,1);
% 
% v0 = zeros(3,1);
% vd0 = [0 ; -g ; 0];
% 
% w1 = R10 * w0 + qd(1)*[0 0 1].';
% wd1 = R10 * wd0 + R10 * cross(w0,qd(1)*[0 0 1].') + qdd(1)*[0 0 1].';
% 
% vd1 = R10 * (cross(wd0,P01) + cross(w0, cross(w0,P01))+vd0);
% vcd1 = cross(wd1,PC1) + cross(w1,cross(w1,PC1)) + vd1;
% 
% F1 = m(1) * vcd1;
% N1 = IC1*wd1 + cross(w1,IC1*w1);
% 
% w2 = R21 * w1 + qd(2)*[0 0 1].';
% wd2 = R21 * wd1 + R21 * cross(w1,qd(2)*[0 0 1].') + qdd(2)*[0 0 1].';
% 
% vd2 = R21 * (cross(wd1,P12) + cross(w1,cross(w1,P12)) + vd1);
% vcd2 = cross(wd2,PC2) + cross(w2, cross(w2,PC2)) + vd2;
% 
% F2 = m(2)*vcd2;
% N2 = IC2*wd2 + cross(w2,IC2*w2);
% 
% f2 = R23 * f3 + F2;
% n2 = N2 + R23 * n3 + cross(PC2,F2) + cross(P23, R23*f3);
% 
% f1 = R12 * f2 + F1;
% n1 = N1 + R12*n2 + cross(PC1,F1) + cross(P12,R12*f2);
% tau1 = transpose(n1) * [0 0 1]';
% tau2 = transpose(n2) * [0 0 1]';
% TAU  = [tau1;tau2];


%% Helper Functions


% 2R IK

function q = IK(x,L,flag)
    Xe = x(1); Ye = x(2); 
    L1 = L(1); L2 = L(2);
    % X/Y Positioning
    length = sqrt(Xe^2+Ye^2);
    if length^2>(L1+L2)^2||length^2<(L2-L1)^2
        theta1 = NaN;
        theta2 = NaN;
        % fprintf('\n Not in Workspace X/Y \n')
        q = [theta1;theta2];
        return
    end
    theta2_up = real(pi-acos((L1^2 + L2^2 - length^2)/(2*L1*L2)));
    theta2_down = real(-(pi-acos((L1^2 + L2^2 - length^2)/(2*L1*L2))));
    theta1_up = atan2(Ye,Xe) - atan2(L2*sin(theta2_up),L1+L2*cos(theta2_up));
    theta1_down = atan2(Ye,Xe) - atan2(L2*sin(theta2_down),L1+L2*cos(theta2_down));
    if flag == 1 
        theta1 = theta1_up;
        theta2 = theta2_up;
        % fprintf('\n Successful Elbow Up! \n')
    else
        theta1 = theta1_down;
        theta2 = theta2_down;
        % fprintf('\n Successful Elbow Down! \n')
    end    
    q = [theta1;theta2];
end

function J = fastJ(q,L)
    J = [ -L(1)*sin(q(1)) - L(2)*sin(q(1)+q(2)) , -L(2)*sin(q(1)+q(2));
           L(1)*cos(q(1)) + L(2)*cos(q(1)+q(2)) ,  L(2)*cos(q(1)+q(2));
           1                                    ,  1                  ;];
end

function Jv = fastJv(q,L)
    Jv = [ -L(1)*sin(q(1)) - L(2)*sin(q(1)+q(2)) , -L(2)*sin(q(1)+q(2));
           L(1)*cos(q(1)) + L(2)*cos(q(1)+q(2)) ,  L(2)*cos(q(1)+q(2));];
end

function Jv_d = fastJv_dot(q,qd,L)
    Jv_d = [ -L(1)*cos(q(1))*qd(1) - L(2)*cos(q(1)+q(2))*(qd(1)+qd(2)) , -L(2)*cos(q(1)+q(2))*(qd(1)+qd(2));
          -L(1)*sin(q(1))*qd(1) - L(2)*sin(q(1)+q(2))*(qd(1)+qd(2)) , -L(2)*sin(q(1)+q(2))*(qd(1)+qd(2));];
end

function P = fastP(q,L)
    P = [L(2)*cos(q(1) + q(2)) + L(1)*cos(q(1));
         L(1)*sin(q(1)) + L(2)*sin(q(1) + q(2))];
end

function q = fastIK(x,L)
    Xe = x(1); Ye = x(2); 
    L1 = L(1); L2 = L(2);
    % X/Y Positioning
    length = sqrt(Xe^2+Ye^2);
    theta2 = real(-(pi-acos((L1^2 + L2^2 - length^2)/(2*L1*L2))));
    theta1= atan2(Ye,Xe) - atan2(L2*sin(theta2),L1+L2*cos(theta2));
    q = [theta1;theta2];
end

% DH Matrix
function T = DH_Matrix(dh)
    alpha = dh(1);
    theta = dh(2);
    a     = dh(3);
    d     = dh(4); 
    T = [ cos(theta)            , -sin(theta)            ,  0          ,  a             ;
          sin(theta)*cos(alpha) ,  cos(theta)*cos(alpha) , -sin(alpha) , -d * sin(alpha);
          sin(theta)*sin(alpha) ,  cos(theta)*sin(alpha) ,  cos(alpha) ,  d * cos(alpha);
          0                     ,  0                     ,  0          ,  1             ;];
end

% Euler Angles from Rotation Matrix
function angles = symT2EUL(T)
    rot = T(1:3,1:3);
    theta1 = simplify(-sin(rot(3,1)),'Steps',15);
    theta2 = simplify(pi-theta1);
    if simplify(rot(3,1)==1) || simplify(rot(3,1)==-1)
        phi = 0;
        if simplify(rot(3,1)==1)
            theta = pi/2;
            psi = simplify(phi + atan(rot(1,2)/rot(1,3)),'Steps',15);
        else 
            theta = -pi/2;
            psi = simplify(-phi + atan(-rot(1,2)/-rot(1,3)),'Steps',15);
        end
        angles = simplify([theta;psi;phi],'Steps',15);
    else
        psi1 = simplify(atan((rot(3,2)/cos(theta1))/(rot(3,3)/cos(theta1))),'Steps',15); 
        psi2 = simplify(atan((rot(3,2)/cos(theta2))/(rot(3,3)/cos(theta2))),'Steps',15); 
        phi1 = simplify(atan((rot(2,1)/cos(theta1))/(rot(1,1)/cos(theta1))),'Steps',15);
        phi2 = simplify(atan((rot(2,1)/cos(theta2))/(rot(1,1)/cos(theta2))),'Steps',15);
        angles = simplify([theta1;psi1;phi1],'Steps',15);
    end
end
March 18 at 6:24 PM
