% loading the dataset
data = load("data.mat").data;

% x,y,z=0, alpha=0, gamma=0, theta
posesg = data.posesg;
poses = data.poses;

new_posesg = [posesg(1:2, :); zeros(1, size(posesg, 2));  zeros(2, size(posesg, 2)); posesg(3, :)];
new_poses = [poses(1:2, :); zeros(1, size(poses, 2));zeros(2, size(poses, 2)); poses(3, :)];

%  x,y,theta=0
landmarksg = data.landmarksg;
landmarks = data.landmarks;

new_landmarksg = [landmarksg; zeros(1, size(landmarksg, 2))];
new_landmarks = [landmarks; zeros(1, size(landmarks, 2))];

% Updating measurements of observations
for i = 1:length(data.observations)
    for j = 1:length(data.observations(i).observation)
        % current measurement
        current_measurement = data.observations(i).observation(j).measurement;
     
        new_measurement = [current_measurement; 0];
     
        % Updating the measurement 
        data.observations(i).observation(j).measurement = new_measurement;
    end
end


% Updating transitions
for i= 1:length(data.transitions)
    v = data.transitions(i).v;
    new_v = [v(1:2, :); 0; 0;0; v(3, :)];
    data.transitions(i).v = new_v;
end

data.posesg = new_posesg;
data.poses = new_poses;
data.landmarksg = new_landmarksg;
data.landmarks = new_landmarks;
 
save('3d_data.mat', 'data');