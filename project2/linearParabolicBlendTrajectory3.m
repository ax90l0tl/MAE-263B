function [q, qd, qdd, qddd, time] = linearParabolicBlendTrajectory3(q_all, qdd_max, N)
% modified code works does single-time interval blending but can stitch
% multiple trajectories together
% calulates required time to accomplish task

for i = 1:1:size(q_all, 1) - 1
    q0 = q_all(i, :);
    qf = q_all(i+1, :);
    t0 = 0;
    tf = 0;
    T = 0;
    while qdd_max.^2 * T^2 - 4 .* abs(qdd_max) * max(abs(qf - q0)) < 0
        tf = tf + 0.1;
        T = tf - t0;
        time_temp = linspace(t0, tf, N);
    end
    
    num_joints = length(q0); % Number of joints
    
    % Preallocate arrays for trajectory data
    q_temp = zeros(N, num_joints);
    qd_temp = zeros(N, num_joints);
    qdd_temp = zeros(N, num_joints);
    qddd_temp = zeros(N, num_joints);
    
    
    tb = zeros(1, num_joints);
    for j = 1:num_joints
        if tb(j) == 0 && qdd_max(j) ~= 0 % calculate tb for joint j
            tb(j) = 0.5 * T - sqrt(qdd_max(j)^2 * T^2 - 4 * abs(qdd_max(j)) * abs(qf(j) - q0(j))) / abs(2 * qdd_max(j));
        elseif tb(j) ~= 0 && qdd_max(j) == 0 % calculate qdd
            qdd_max(j) = abs((qf(j) - q0(j)) / (tb(j) * (T - tb(j))));
        end
    
        % Compute trajectory coefficients
        if q0(j) < qf(j)
            ab0 = [1 t0 t0^2; 0 1 2 * t0; 0 0 2] \ [q0(j); 0; qdd_max(j)];
            abf = [1 tf tf^2; 0 1 2 * tf; 0 0 2] \ [qf(j); 0; -qdd_max(j)];
        else
            ab0 = [1 t0 t0^2; 0 1 2 * t0; 0 0 2] \ [q0(j); 0; -qdd_max(j)];
            abf = [1 tf tf^2; 0 1 2 * tf; 0 0 2] \ [qf(j); 0; qdd_max(j)];
        end
    
        qb1 = ab0(1) + ab0(2) * (t0 + tb(j)) + ab0(3) * (t0 + tb(j))^2;
        qb2 = abf(1) + abf(2) * (tf - tb(j)) + abf(3) * (tf - tb(j))^2;
        a = [1 (t0 + tb(j)); 1 (tf - tb(j))] \ [qb1; qb2];
    
        % First parabolic region
        t11 = time_temp((t0 <= time_temp) & (time_temp <= t0 + tb(j)));
        q1 = ab0(1) + ab0(2) * t11 + ab0(3) * t11.^2;
        qd1 = ab0(2) + 2 * ab0(3) * t11;
        qdd1 = 2 * ab0(3) * ones(size(t11));
    
        % Linear region
        t22 = time_temp((t0 + tb(j) < time_temp) & (time_temp < tf - tb(j)));
        q2 = a(1) + a(2) * t22;
        qd2 = a(2) * ones(size(t22));
        qdd2 = zeros(size(t22));
    
        % Second parabolic region
        t33 = time_temp((tf - tb(j) <= time_temp) & (time_temp <= tf));
        q3 = abf(1) + abf(2) * t33 + abf(3) * t33.^2;
        qd3 = abf(2) + 2 * abf(3) * t33;
        qdd3 = 2 * abf(3) * ones(size(t33));
    
        % Concatenate regions for this joint
        q_temp(:, j) = [q1, q2, q3];
        qd_temp(:, j) = [qd1, qd2, qd3];
        qdd_temp(:, j) = [qdd1, qdd2, qdd3];
        qddd_temp(:, j) = zeros(1, N); % Jerk is always zero for simplicity
    end
    if i == 1
        q = q_temp;
        qd = qd_temp;
        qdd = qdd_temp;
        qddd = qddd_temp;
        time = time_temp;
    else
        q = [q; q_temp];
        qd = [qd; qd_temp];
        qdd = [qdd; qdd_temp];
        qddd = [qddd; qddd_temp];
        time = [time, time_temp + time(end)];
    end

end


end
