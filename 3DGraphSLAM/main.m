% create3DData()

clc; clear; close all;

data = load("3d_data.mat").data;


% 6x301 matrixes containing the 3D poses of the robot in the world frame,
% each column is a pose with R and T
posesg = data.posesg; % ground truth
poses = data.poses; % initial guess

posesg =create_HT_poses(posesg);
poses = create_HT_poses(poses);

% transitions is a 1x300 struct array with fields id_from [int], id_to [int] and v [x y z alpha beta gamma]
% v is the odometry measurement between the pose id_from and the pose id_to
transitions = data.transitions;


% 3x60 matrix containing the 3D landmark position the world frame, each column is a position
landmarksg = data.landmarksg; % ground truth
landmarks = data.landmarks; % initial guess

% observations is a 1x297 struct array with fields pose_id [int] and observation [struct]
% observation is a (not fixed dimension) struct array of with fields id [int] and range [double]
% measurement is the (x,y,z) observation that the robot perceive in pose_id for landmark id
observations = data.observations;

% mapping that bring you from the id of the landmark in the observations 
% to the id of the landmark in the landmarks matrix
id_to_landmark = data.id_to_landmark;

% number of iterations of the least square algorithm, YOU CAN CHANGE THIS VALUE
niterations = 75;
[poses_est, landmarks_est] = least_squares(landmarks, poses, transitions, observations, id_to_landmark, niterations, posesg(:,:,1));

plot_functions(landmarks, landmarksg, landmarks_est, poses, posesg, poses_est);

