function [q,qd, qdd, qddd, time] = linearParabolicBlendTrajectory2(traj, t, v_max, a_max, N)
% My attempt at doing multi-via point blending
% Works for some cases but really particular about timing and its easy to
% break
% traj is a n x ndof array of joint angles along a trajectory
% v_max is max velocity for each joint
% a_max is max acceleration for each joint
% N is density of steps between each set of points

ndof = size(traj, 2);
num_pts = size(traj, 1);
min_n = 50; % minimum density of trajectory
t_blend = t;
qd_jk = 0;
% calculate blending times
for i = 1:1:num_pts
    if i == 1
        q_0 = traj(i, :);
        q_1 = traj(i + 1, :);
        dq = q_1 - q_0;
        qdd = ((dq == 0) * 1 + (dq ~= 0) .* sign(dq)) .* a_max;
        dt = t(i);
        t1 = dt - sqrt(dt.^2 - 2.*(dq)./qdd);
        assert(all(t1 >= 0))
        t_blend(i, :) = t1;
        qd_jk = dq./(dt - 0.5.*t1);
    else
        if i == num_pts
            q_0 = traj(i - 1, :);
            q_1 = traj(i, :);
            dt = t(i-1);
            dq = q_0 - q_1;
            qdd = ((dq == 0) * 1 + (dq ~= 0) .* sign(dq)) .* a_max;
            tn = dt - sqrt(dt.^2 - 2.*(dq)./qdd);
            t_blend(i, :) = tn;
        end
        if i ~= 1 && i ~= num_pts
            q_0 = traj(i, :);
            q_1 = traj(i+1, :);
            dq = q_1 - q_0;
            dt = t(i);
            qd_kl = dq./dt;
            qdd = ((dq == 0) * 1 + (dq ~= 0) .* sign(dq)) .* a_max;
            tk = (qd_kl - qd_jk)./(sign(qd_kl - qd_jk).*a_max);
            t_blend(i, :) = tk;
            qd_jk = qd_kl;
        end
    end
end

% calculate linear portion times
t_line = zeros(size(t_blend, 1) - 1, ndof);
for i = 1:1:size(t_blend, 1) - 1
    if i == 1
        t_line(i, :) =  t(i, :) - t_blend(i, :) - 0.5*t_blend(i + 1, :);
    elseif i == size(t_blend, 1) - 1
        t_line(i, :) =  t(i, :) - t_blend(i+1, :) - 0.5*t_blend(i, :);
    else
        t_line(i, :) =  t(i, :) - 0.5*t_blend(i, :) - 0.5*t_blend(i + 1, :);
    end
end
% t_blend
% t_line
% calculate q, qd, qdd, qddd
for i = 1:1:num_pts
    % last parabolic blend
    % ensures velocity ends at 0
    if i == num_pts
        q_0 = traj(i-1, :);
        q_1 = traj(i, :);
        dq = q_1 - q_0;
        qd_temp = qd(end, :);
        dt_i = t_blend(i, :);
        n = max([round((dt_i./dt).*N), min_n]);
        t_i = transpose(linspace(0, 1, n)) * dt_i;
        non_zero_indices = t_i(end, :) ~= 0;
        qdd_temp = zeros(1, ndof);
        qdd_temp(non_zero_indices) = -qd_temp(non_zero_indices) ./ t_i(end, non_zero_indices);
        % assert(all(abs(qdd_temp) <= a_max))
        qdd_temp = min(a_max, qdd_temp);
        q_temp = q(end, :) + qd_temp.*(t_i) + 0.5*qdd_temp.*(t_i).^2;
        q = [q; q_temp];
        qd = [qd; qd_temp + qdd_temp.*(t_i)];
        qdd = [qdd; qdd_temp .* ones(n, ndof)];
        qddd = [qddd; zeros(n, ndof)];
        time = [time; t_i + time(end, :)];

        % do one last parabolic blend to ensure we get to our destination
        % dont have time to make this function handle multi arrays
        
        q_f = zeros(min_n, ndof);
        qd_f = zeros(min_n, ndof);
        qdd_f = zeros(min_n, ndof);
        qddd_f = zeros(min_n, ndof);
        time_f = zeros(min_n, ndof);
        for j = 1:1:ndof
            [qf_t, qdf_t, qddf_t, qdddf_t, timef_t] = linearParabolicBlendTrajectory(time(end, j), time(end, j), q(end, j), q_1(j), 0, a_max(j), min_n);
            q_f(:, j) = transpose(qf_t);
            qd_f(:, j) = transpose(qdf_t);
            qdd_f(:, j) = transpose(qddf_t);
            qddd_f(:, j) = transpose(qdddf_t);
            time_f(:, j) = transpose(timef_t);
        end
        q = [q; q_f];
        qd = [qd; qd_f];
        qdd = [qdd; qdd_f];
        qddd = [qddd; qddd_f];
        time = [time; time_f];
        break;
    end
    q_0 = traj(i, :);
    q_1 = traj(i + 1, :);
    dq = q_1 - q_0;
    qdd_temp = sign(dq).*a_max;
    dt = t(i, :);
    dt_i = t_blend(i, :);
    dt_line = t_line(i, :);
    
    % first parabolic blend
    if i == 1
        n = max([round((dt_i./dt).*N), min_n]);
        t_i = transpose(linspace(0, 1, n)) * dt_i;
        q_para_1 = q_0 + 0.5*qdd_temp.*t_i.^2;
        qd_temp = qdd_temp.*t_i;

        q = q_para_1;
        qd = qd_temp;
        qdd = qdd_temp .* ones(n, ndof);
        qddd = zeros(n, ndof);

        time = t_i;

    end
    % Linear sequences
    n = max([round((dt_line./dt).*N), min_n]);
    t_lin = transpose(linspace(0, 1, n)) * dt_line;
    if i == 1
        qd_temp = dq./(dt-0.5.*dt_i);
    else
        qd_temp = dq/dt;
    end
    q = [q; q(end, :) + qd_temp.*(t_lin)];
    qd = [qd; qd_temp.*ones(n, ndof)];
    qdd = [qdd; zeros(n, ndof)];
    qddd = [qddd; zeros(n, ndof)];
    time = [time; t_lin + time(end, :)];
    
    if i+1 == num_pts
        continue
    end

    % Parabolic blending
    dt_i = t_blend(i+1, :);
    n = max([round((dt_i./dt).*N), min_n]);
    t_i = transpose(linspace(0, 1, n)) * dt_i;

    qdd_temp = sign((traj(i + 2, :) - traj(i+1, :)./t(i+1, :)) - qd_temp).*a_max;

    q_temp = q(end, :) + qd_temp.*(t_i) +0.5.*qdd_temp.*(t_i).^2;
    qd_temp = qd_temp + qdd_temp.*(t_i);
    q = [q; q_temp];
    qd = [qd; qd_temp];
    qdd = [qdd; qdd_temp .* ones(n, ndof)];
    qddd = [qddd; zeros(n, ndof)];
    time = [time; t_i + time(end, :)];

end




% optimize shortest time blend (trianguler velocity)
% figure();
% hold on
% qd_0 = 0;
% q_0 = traj(end -1);
% q = traj(end);
% qdd = sign(q - q_0)*a_max;
% t1 = linspace(0, sqrt(abs(q)/a_max), N/2);
% t2 = linspace(sqrt(abs(q)/a_max), 2*sqrt(abs(q)/a_max), N/2);
% A_q = q_0 + qd_0*t1 + 0.5*qdd*t1.^2;
% A_qd = qd_0 + qdd*t1;
% B_q = A_q(end) + A_qd(end)*(t2 - t1(end)) - 0.5*qdd*(t2-t1(end)).^2;
% B_qd = A_qd(end) - qdd*(t2 - t1(end));
%
% plot(t1, A_q)
% plot(t2, B_q)
%
% figure();
% hold on
% plot(t1, A_qd)
% plot(t2, B_qd)
% t1_blend = linspace(0, t1, )

end