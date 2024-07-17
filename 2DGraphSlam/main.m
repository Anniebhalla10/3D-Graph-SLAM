% dont changeng just change the number of iterations

clc; clear; close all;
% loading the dataset
data = load("data.mat").data;

% 3x301 matrixes containing the 2D poses of the robot in the world frame, each column is a pose
posesg = data.posesg; % ground truth
poses = data.poses; % initial guess
% 2x60 matrix containing the 2D landmark position the world frame, each column is a position
landmarksg = data.landmarksg; % ground truth
landmarks = data.landmarks; % initial guess
% transitions is a 1x300 struct array with fields id_from [int], id_to [int] and v [x y theta]
%   v is the odometry measurement between the pose id_from and the pose id_to
transitions = data.transitions;
% observations is a 1x297 struct array with fields pose_id [int] and observation [struct]
%   observation is a (not fixed dimension) struct array of with fields id [int] and range [double]
%   measurement is the (x,y) observation that the robot perceive in pose_id for landmark id
observations = data.observations;
% mapping that bring you from the id of the landmark in the observations 
% to the id of the landmark in the landmarks matrix
id_to_landmark = data.id_to_landmark;

global all_poses all_landmarks

all_poses={};
all_landmarks={};

all_poses{1}= poses;
all_landmarks{1}= landmarks;

% number of iterations of the least square algorithm, YOU CAN CHANGE THIS VALUE
niterations = 23;
[poses_est, landmarks_est] = least_square(landmarks, poses, transitions, observations, id_to_landmark, niterations, posesg(:,1:3));
% plot of the trajectories
plot_trajectories(poses', poses_est', posesg')
saveas(gcf,'Trajectories.png')
% plot of the maps (landmarks)
plot_maps(landmarks', landmarks_est', landmarksg', false)
saveas(gcf,'Landmarks.png')
% both the trajectories and the maps are plotted in the same figure
plot_trajectories(poses', poses_est', posesg')
plot_maps(landmarks', landmarks_est', landmarksg', true)
saveas(gcf,'TrajectoriesAndLandmarks.png')

plot_all(all_poses, all_landmarks)

function plot_maps(landmarks, landmarksest, landmarksg, holdon)
    if holdon, hold on; title("Full result"); else, figure(); hold on; axis equal; title("Maps comparison"); end
    ylim([-15 12]);
    plot(reshape(landmarks(:, 1), 1, []), reshape(landmarks(:, 2), 1, []), 'o','MarkerEdgeColor','red', 'MarkerSize', 10,'LineWidth',2);   
    plot(reshape(landmarksest(:, 1), 1, []), reshape(landmarksest(:, 2), 1, []), 'o','MarkerEdgeColor','blue', 'MarkerSize', 10,'LineWidth',4);   
    plot(reshape(landmarksg(:, 1), 1, []), reshape(landmarksg(:, 2), 1, []), 'o','MarkerEdgeColor','green', 'MarkerSize', 10,'LineWidth',2);
    legend("initial guess", "after ls", "ground truth");
end

function plot_trajectories(poses, posesest, posesg)
    figure(); axis equal; hold on;
    title("Trajectories comparison"); 
    plot(reshape(posesest(:, 1), 1, []), reshape(posesest(:, 2), 1, []), 'b-', 'linewidth', 5);   
    plot(reshape(posesg(:, 1), 1, []), reshape(posesg(:, 2), 1, []), 'g-', 'linewidth', 3);   
    plot(reshape(poses(:, 1), 1, []), reshape(poses(:, 2), 1, []), 'r-', 'linewidth', 3);  
    legend("initial guess", "after ls", "ground truth");
end

function plot_all(poses, landmarks)
    figure();
    hold on;
    axis equal;
    axis([-15 15 -15 15]);
    for i = 1:length(poses)
            title("FUll evolution: iteration "+(i-1));
        if i~=1
            delete(h1)
            delete(h2)
            axis([-15 15 -15 15]);
        end
        po= poses{i}';
        la= landmarks{i}';
        h1= plot(reshape(po(:,1),1 ,[]), reshape(po(:,2),1,[]), 'c-', 'LineWidth',5);

        h2= plot(reshape(la(:,1),1 ,[]), reshape(la(:,2),1,[]), 'o', 'MarkerEdgeColor','magenta','MarkerSize', 10,'LineWidth',2);
        pause(0.5);
    end

end
