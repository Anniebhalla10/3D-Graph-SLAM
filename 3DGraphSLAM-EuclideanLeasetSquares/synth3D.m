%% Create a synthetic 3D dataset for PLGO (Pose Landmark Graph Optimization)

% 0. Dataset parameters
n_poses = 20;
n_landmarks = 100;

% plot the dataset?
viz = true;

% seed
rng(102030)

R = 10;
syms t real

% 1. Generate a set of poses

% we use a spherical trajectory

% define the azimuth and elevation trajectory

% (singularities)
az = 2*sin(3*t);

% no singularities
% az = 2*sin(3*t+ 0.3);
el = t;

% and map it to Cartesian
x = R * cos(el) .* cos(az);
y = R * cos(el) .* sin(az);
z = R * sin(el);

t_ = linspace(0,2*pi, n_poses);
axisLength = 2;

% Poses as 6D Euclidean Vector
% HT as homogeneous transformations
HT = zeros(4,4,n_poses);
posesg= zeros(6, n_poses);
for i = 1:length(t_)
    
    % Translation part
    pos = eval(subs([x;y;z], t, t_(i)));
    
    % Rotation part
    tangent = eval(subs(jacobian([x;y;z], t), t_(i)));
    xAxis = tangent / norm(tangent);
    
    zAxis = pos / norm(pos);
    
    yAxis = cross(zAxis, xAxis);
    yAxis = yAxis / norm(yAxis);
    
    xAxis = cross(yAxis, zAxis);
    
    HT(:,:,i)=[xAxis yAxis zAxis pos;
        0     0     0   1];
    
    posesg(:,i) = extract6DPose(HT(:,:,i));
end

% Remove singular poses
% [posesg, HT] = pruneSingularityPoses(posesg, HT);

% function [posesg, HT] = pruneSingularityPoses(posesg, HT)
%     disp("here")
%     indicesToRemove = all(posesg == 0, 1);
%     posesg(:, indicesToRemove) = [];
%     HT(:, :, indicesToRemove) = [];
% end

% 2. Generate a set of landmarks ground truth

% we once again make use of spherical coordinates
% intervals for the points;
r_int = [R+2 R+8];
az_int = [0 2*pi];
el_int = [0 2*pi];

limits = [r_int; az_int; el_int];

l_sph = (limits(:,2) - limits(:,1)).*rand(3,100) + limits(:,1) ;

landmarksg = [l_sph(1,:) .* cos(l_sph(3,:)) .* cos(l_sph(2,:));
    l_sph(1,:) .* cos(l_sph(3,:)) .* sin(l_sph(2,:));
    l_sph(1,:) .* sin(l_sph(3,:)) ] ;


% 3. generate pose-pose perfect measurements
transition = struct("from", -1, "to", -1, "meas", zeros(6,1));
transitions = transition;
for i=1:size(posesg, 2)-1
    transitions(i) = transition;
    transitions(i).from = i;
    transitions(i).to = i+1;
    z_ = HT(:,:,i)\HT(:,:,i+1);
    transitions(i).meas = extract6DPose(z_);
end

% 4. generate pose-landmark perfect measurements
observation = struct("p_idx", -1, "l_idx", -1, "meas", zeros(3,1));
observations = observation([]);
% we iterate over landmarks and assign to them N poses that see this
% landmark
N = 3;
for i = 1:size(landmarksg,2)
    poses_idxes = randperm(size(posesg,2), N);  % Random pose indices
    for j = 1:N
        observations(end+1).p_idx = poses_idxes(j);
        observations(end).l_idx = i;
        rel_measurement = HT(:,:,poses_idxes(j)) \ [landmarksg(:,i); 1];
        observations(end).meas = rel_measurement(1:3);
    end
end

% 5. add noise to the measurements
% Add noise to pose-pose measurements (transitions)
pose_noise_std = 0.01; % standard deviation for pose noise
for i = 1:length(transitions)
    % Extract current poses from transition
    from_pose = transitions(i).from;
    to_pose = transitions(i).to;
    
    noise_trans = pose_noise_std * randn(3, 1);
    noise_rot = pose_noise_std * randn(3, 1);
    
    % Apply noise to the measurement
    transitions(i).meas(1:3) = transitions(i).meas(1:3) + noise_trans;
    transitions(i).meas(4:6) = transitions(i).meas(4:6) + noise_rot;
end

% Add noise to pose-landmark measurements (observations)
landmark_noise_std = 0.01; % standard deviation for landmark noise
for i = 1:length(observations)
    noise = landmark_noise_std * randn(3, 1);
    observations(i).meas = observations(i).meas + noise;
end

% 6. Generate intial guesses for the landmarks and poses
landmark_noise_std = 0.9;
landmarks = landmarksg + landmark_noise_std * randn(size(landmarksg));

pose_noise_std = 0.9;
poses = zeros(6, size(posesg,2));
for i = 1:size(posesg,2)
    current_pose = posesg(:, i);
    
    % Add noise to the translational part
    noise_trans = pose_noise_std * randn(3, 1);
    
    % Add noise to the rotational part
    noise_rot = pose_noise_std * randn(3, 1);
    
    % Apply noise to position and rotation parts of the 6D pose
    poses(1:3, i) = current_pose(1:3) + noise_trans;
    poses(4:6, i) = current_pose(4:6) + noise_rot;
end

% save the dataset
save("synth3D.mat", "posesg", "landmarksg", "transitions", "observations", "poses", "landmarks")

%plot the dataset
if viz
    [xs, ys, zs] = sphere;
    xs = R * xs;
    ys = R * ys;
    zs = R * zs;
    figure;
    mesh(xs, ys, zs, 'EdgeColor', [0.5, 0.5, 0.5]); hold on; axis equal;
    
    fplot3(x, y, z, 'b', 'LineWidth', 2);
    for i=1:size(posesg,2)
        hold on;
        plotPose(posesg(:,i))
    end
    
    scatter3(landmarksg(1,:), landmarksg(2,:), landmarksg(3,:));
end

function hnd = plotPose(pose)
hold on;

pos = pose(1:3);
R = Rx(pose(4))*Ry(pose(5))*Rz(pose(6));

axisLength = 2;

quiver3(pos(1), pos(2), pos(3), R(1,1)*axisLength, R(2,1)*axisLength, R(3,1)*axisLength, 'color', 'r', 'LineWidth', 2);
quiver3(pos(1), pos(2), pos(3), R(1,2)*axisLength, R(2,2)*axisLength, R(3,2)*axisLength, 'color', 'g', 'LineWidth', 2);
quiver3(pos(1), pos(2), pos(3), R(1,3)*axisLength, R(2,3)*axisLength, R(3,3)*axisLength, 'color', 'b', 'LineWidth', 2);
end