function [q,qd,qdd, qddd] = linearParabolicBlendTrajectory_inclass(t0, tf, q0, qf, tb, qdd, N)% linearParabolicBlendTrajectory:
% linearParabolicBlend generates joint-space positions from q0 to qf
%   with a trapezoidal velocity profile (accel-constant-accel).
%
%   INPUTS:
%       q0  - initial joint angles 
%       qf  - final joint angles   
%       t0  - initial time 
%       tf  - final time
%       tb  - blend time (calculate qd)
%       qdd - desired acceleration (calculate blend time)

%       N   - number of trajectory points
%
% Outputs:
%   pos - [1 x N] positions
%   vel - [1 x N] velocities
%   acc - [1 x N] accelerations
%   jerk   - [1 x N] jerk

T = tf - t0;
% Time vector
t = linspace(t0, tf, N);

if(tb == 0 && qdd ~= 0) % calculate tb
    tb = [];
elseif(tb ~= 0 && qdd == 0) % calculate qdd
    qdd = [];
end

% if qdd < 4 * (qf - q0) / T^2 return  !!!

    
if q0 < qf
    ab0 = [] \ [q0 0 qdd]'; % Coefficient for first blend
    abf = [] \ [qf 0 -qdd]'; % Coefficient for second blend
else
    ab0 = [] \ []; % Coefficient for first blend
    abf = [] \ []; % Coefficient for second blend
end
qb1 = []
qb2 = []

a = [] % Coefficient for linear region
% first parabolic region
t11 = t((t0<=t) & (t<=t0+tb));
q = []
qd  = []
qdd = [];
% Linear region
t22 = t((t0+tb<t) & (t<tf-tb)); % linear region
q = [q ]
qd = [qd ]
qdd = [qdd ]
% second parabolic region
t33 = t((tf-tb<=t) & (t<=tf)); 
q   = [q ]
qd  = [qd ]
qdd = [qdd ]
qddd = zeros(1, N);
end
