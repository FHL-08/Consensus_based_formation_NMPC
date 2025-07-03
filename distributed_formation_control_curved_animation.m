clear
clc
close all

load('Curved_Trajectory_Values.mat');

%%% Define video settings %%%
videoFile = 'distributed_formation_control_curved.avi';  % File name
v = VideoWriter(videoFile, 'Motion JPEG AVI'); % Create video writer
v.FrameRate = 30; % Set frame rate (adjust as needed)
open(v); % Open the video file for writing

robot_radius = 0.2;
virtual_path = animatedline('Color', 'c', 'LineStyle', '--', 'LineWidth', 1.5);
trail1 = animatedline('Color', 'r', 'LineWidth', 1.5);
trail2 = animatedline('Color', 'b', 'LineWidth', 1.5);
trail3 = animatedline('Color', 'g', 'LineWidth', 1.5);

f=figure(1);
hold on;
%%% Obstacle %%%
rectangle('Position', [p_obs(1)-obstacle_radius,...
                   p_obs(2)-obstacle_radius,...
                   2*obstacle_radius, 2*obstacle_radius],...
                  'Curvature', [1, 1], 'FaceColor',[0 0 0]);
obstacle_virtual = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);

%%% Robot 1 %%%
robot1_body = rectangle('Position', [q_i_library(1,1,1)-robot_radius,...
                   q_i_library(2,1,1)-robot_radius,...
                   2*robot_radius, 2*robot_radius],...
                  'Curvature', [1, 1], 'FaceColor',[1 0 0]);
robot1_heading = plot([q_i_library(1,1,1)+robot_radius*cos(q_i_library(3,1,1));q_i_library(1,1,1)],...
                      [q_i_library(2,1,1)+robot_radius*sin(q_i_library(3,1,1));q_i_library(2,1,1)],...
                      'y','LineWidth',1.5);
robot1_virtual = plot(nan, nan, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);

%%% Robot 2 %%%
robot2_body = rectangle('Position', [q_i_library(1,1,2)-robot_radius,...
                   q_i_library(2,1,2)-robot_radius,...
                   2*robot_radius, 2*robot_radius],...
                  'Curvature', [1, 1], 'FaceColor',[0 0 1]);
robot2_heading = plot([q_i_library(1,1,2)+robot_radius*cos(q_i_library(3,1,2));q_i_library(1,1,2)],...
                      [q_i_library(2,1,2)+robot_radius*sin(q_i_library(3,1,2));q_i_library(2,1,2)],...
                      'y','LineWidth',1.5);
robot2_virtual = plot(nan, nan, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10);

%%% Robot 3 %%%
robot3_body = rectangle('Position', [q_i_library(1,1,3)-robot_radius,...
                   q_i_library(2,1,3)-robot_radius,...
                   2*robot_radius, 2*robot_radius],...
                  'Curvature', [1, 1], 'FaceColor',[0 1 0]);
robot3_heading = plot([q_i_library(1,1,3)+robot_radius*cos(q_i_library(3,1,3));q_i_library(1,1,3)],...
                      [q_i_library(2,1,3)+robot_radius*sin(q_i_library(3,1,3));q_i_library(2,1,3)],...
                      'y','LineWidth',1.5);
% plot(q_i_library(1, 1:i, 3), q_i_library(2, 1:i, 3), 'g', 'LineWidth', 1.5);
robot3_virtual = plot(nan, nan, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

%%% Formation Shape %%%
cord1x = reshape(q_i_library(1, 1, 1:3), [], 1);
cord1y = reshape(q_i_library(2, 1, 1:3), [], 1);
formation_shape = plot([cord1x; cord1x(1); p_l(1, 1)], [cord1y; cord1y(1); p_l(2, 1)], '.-',...
                    'Color', [0.5, 0.5, 0.5, 0.5],'MarkerSize', 10, 'LineWidth', 2);

%%% FORMATTING %%%
legend('show', 'Location', 'northeast');
legend([virtual_path, obstacle_virtual, robot1_virtual, robot2_virtual, robot3_virtual],...
    {'Virtual Leader', 'Obstacle', 'Robot 1', 'Robot 2','Robot 3'});
title("Robots' Trajectories", 'fontsize', 15);
xlabel('X coordinate [m]', 'fontSize', 12);
xlim([-3 10]);
ylabel('Y coordinate [m]', 'fontSize', 12);
ylim([-1.5 12]);
pbaspect([diff(xlim) diff(ylim) 1]);
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');
set(f, 'Position', [100, 100, 1200, 800]);

%%% Animation %%%
for i=1:length(k)
    if ~isvalid(f) || ~ishandle(f)
        disp('Figure closed. Animation cancelled.');
        break;
    end

    %%% Update Leader Path %%%
    addpoints(virtual_path, p_l(1, i), p_l(2, i));

    %%% Update Robot 1 %%%
    robot1_body.Position = [q_i_library(1,i,1)-robot_radius,...
                            q_i_library(2,i,1)-robot_radius,...
                            2*robot_radius, 2*robot_radius];
    robot1_heading.XData = [q_i_library(1,i,1)+robot_radius*cos(q_i_library(3,i,1));q_i_library(1,i,1)];
    robot1_heading.YData = [q_i_library(2,i,1)+robot_radius*sin(q_i_library(3,i,1));q_i_library(2,i,1)];
    addpoints(trail1, q_i_library(1, i, 1), q_i_library(2, i, 1));

    %%% Update Robot 2 %%%
    robot2_body.Position = [q_i_library(1,i,2)-robot_radius,...
                            q_i_library(2,i,2)-robot_radius,...
                            2*robot_radius, 2*robot_radius];
    robot2_heading.XData = [q_i_library(1,i,2)+robot_radius*cos(q_i_library(3,i,2));q_i_library(1,i,2)];
    robot2_heading.YData = [q_i_library(2,i,2)+robot_radius*sin(q_i_library(3,i,2));q_i_library(2,i,2)];
    addpoints(trail2, q_i_library(1, i, 2), q_i_library(2, i, 2));

    %%% Update Robot 3 %%%
    robot3_body.Position = [q_i_library(1,i,3)-robot_radius,...
                            q_i_library(2,i,3)-robot_radius,...
                            2*robot_radius, 2*robot_radius];
    robot3_heading.XData = [q_i_library(1,i,3)+robot_radius*cos(q_i_library(3,i,3));q_i_library(1,i,3)];
    robot3_heading.YData = [q_i_library(2,i,3)+robot_radius*sin(q_i_library(3,i,3));q_i_library(2,i,3)];
    addpoints(trail3, q_i_library(1, i, 3), q_i_library(2, i, 3));

    %%% Update Formation Shape %%%
    cord1x = reshape(q_i_library(1, i, 1:3), [], 1);
    cord1y = reshape(q_i_library(2, i, 1:3), [], 1);
    formation_shape.XData = [cord1x; cord1x(1); p_l(1, i)];
    formation_shape.YData = [cord1y; cord1y(1); p_l(2, i)];
    hold off
    drawnow limitrate;

    %%% Capture frame and write to video %%%
    frame = getframe(f);
    writeVideo(v, frame);
end

% Clsoe the video
close(v);