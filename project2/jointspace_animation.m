function jointspace_animation(q, robot, N, view, filename)
v = VideoWriter(filename,'MPEG-4'); 
v.FrameRate = 20;     
open(v);
p = zeros(N, 3);
for i = 1: N
    T = robot.fkine(q(i,:));  
    p(i,:) = transl(T);     
    plot3(p(i,1), p(i,2), p(i,3),'*r');
    hold on;
    robot.plot(q(i,:), 'view', view, 'workspace', [-1, 1, -1, 1, 0, 1]);
    drawnow;                       % ensure the figure updates
    frame = getframe(gcf);         % capture the figure as a frame
    writeVideo(v, frame);          % write this frame to the video
end
pauseFrames = round(v.FrameRate * 0.8);  % 0.8 second worth of frames

% Capture the last frame (so it stays on screen for the "pause")
lastFrame = getframe(gcf);

for j = 1:pauseFrames
    writeVideo(v, lastFrame);
end

%Close the video
close(v);
end