clc; clear; close all;

data = load("synth3D.mat");

global all_poses all_landmarks chi_values

all_poses={};
all_landmarks={};

posesg = data.posesg;  %ground truth
poses = data.poses;  %initial guess

%transition is a struct containing from (pose id) to (pose id) meas (6D)
transitions = data.transitions;

landmarksg = data.landmarksg;  %ground truth
landmarks = data.landmarks;  %initial guess

%observation struct contains p_idx l_idx and meas (x y z)
observations = data.observations;


all_poses{1}= poses;
all_landmarks{1}= landmarks;

%number of iterations of the least square algorithm, YOU CAN CHANGE THIS VALUE
niterations = 55;

chi_values = zeros(niterations);

[poses_est, landmarks_est] = least_squares(landmarks, poses, transitions, observations, niterations, posesg(:,1));

plotFunction(landmarksg, landmarks, landmarks_est, posesg, poses, poses_est);
