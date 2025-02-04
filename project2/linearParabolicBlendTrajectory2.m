function [q,qd,qdd, qddd] = linearParabolicBlendTrajectory2(q0, qd_max, qdd_max, N)
% linearParabolicBlend generates joint-space positions from q0 to qf
% with a trapezoidal velocity profile (accel-constant-accel).

% INPUTS:
% q0 - initial joint angles
% qf - final joint angles
% t0 - initial time
% tf - final time
% tb - blend time (calculate qd)
% qdd - desired acceleration (calculate blend time)
% N - number of trajectory points
%
% Outputs:
% pos - [1 x N] positions
% vel - [1 x N] velocities
% acc - [1 x N] accelerations
% jerk - [1 x N] jerk

q = zeros(N, size(q0, 2));
qd = zeros(N, size(q0, 2));
qdd = zeros(N, size(q0, 2));
qddd = zeros(N, size(q0, 2));
for i = 1:1:length(q)
    if i == 1
    
    elseif i == length(q)

    else

end

end