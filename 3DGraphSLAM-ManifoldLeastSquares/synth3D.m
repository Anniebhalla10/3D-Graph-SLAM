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
az = 2*sin(3*t);
el = t;

% and map it to Cartesian
x = R * cos(el) .* cos(az);
y = R * cos(el) .* sin(az);
z = R * sin(el);

t_ = linspace(0,2*pi, n_poses);
axisLength = 2;

posesg= zeros(4,4, n_poses);
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
    
    posesg(:,:,i) =[xAxis yAxis zAxis pos; 
                        0     0     0   1];
end


% 2. Generate a set of landmarks groudn truth

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
transition = struct("from", -1, "to", -1, "meas", eye(4));
transitions = transition;
for i=1:size(posesg, 3)-1
    transitions(i) = transition;
    transitions(i).from = i;
    transitions(i).to = i+1;
    transitions(i).meas = posesg(:,:,i)\posesg(:,:,i+1);
end

% 4. generate pose-landmark perfect measurements

observation = struct("p_idx", -1, "l_idx", -1, "meas", zeros(3,1));
observations = observation([]);
% we iterate over landmarks and assign to them N poses that see this
% landmark
N = 3;
for i=1:size(landmarksg,2)
    % the poses that see this landmark
    poses_idxes = randperm(size(posesg,3), N);
    for j=1:N
        observations(end+1).p_idx = poses_idxes(j);
        observations(end).l_idx = i;

        meas_h = posesg(:,:,poses_idxes(j)) \ [landmarksg(:,i); 1];
        observations(end).meas = meas_h(1:3);
    end
end

% 5. add noise to the measurements

% a) Add noise to pose-pose measurements (transitions)
pose_noise_std = 0.01; % standard deviation for pose noise
for i = 1:length(transitions)
    noise_trans = pose_noise_std * randn(3, 1);
    noise_rot = expm(skew(pose_noise_std * randn(3,1)));
    measurement= transitions(i).meas;
    measurement(1:3, 1:3) = measurement(1:3, 1:3) * noise_rot;
    measurement(1:3, 4) = measurement(1:3, 4) + noise_trans;
    transitions(i).meas = measurement;
end

% b) Add noise to pose-landmark measurements (observations)
landmark_noise_std = 0.01; % standard deviation for landmark noise
for i = 1:length(observations)
    noise = landmark_noise_std * randn(3, 1);
    observations(i).meas = observations(i).meas + noise;
end

% 6. Generate intial guesses for the landmarks and poses
landmark_noise_std = 0.9;  
landmarks = landmarksg + landmark_noise_std * randn(size(landmarksg));

pose_noise_std = 0.9; 
poses = zeros(4,4,n_poses);
for i = 1:n_poses
    noise_trans = pose_noise_std * randn(3, 1);
    noise_rot = expm(skew(pose_noise_std * randn(3,1)));
    poses(:,:,i) = posesg(:,:,i);
    poses(1:3, 1:3, i) = posesg(1:3, 1:3, i) * noise_rot;
    poses(1:3, 4, i) = posesg(1:3, 4, i) + noise_trans;
end


% save the dataset
save("synth3D.mat", "posesg", "landmarksg", "transitions", "observations", "poses", "landmarks")

% plot the dataset
if viz
    [xs, ys, zs] = sphere;
    xs = R * xs;
    ys = R * ys;
    zs = R * zs;
    figure;
    mesh(xs, ys, zs, 'EdgeColor', [0.5, 0.5, 0.5]); hold on; axis equal;
    fplot3(x, y, z, 'b', 'LineWidth', 2);
    for i=1:size(posesg,3)
        hold on;
        plotPose(posesg(:,:,i))
    end
    
    scatter3(landmarksg(1,:), landmarksg(2,:), landmarksg(3,:));
end

function hnd = plotPose(pose)
    hold on;

    pos = pose(1:3, 4);
    R = pose(1:3, 1:3);

    axisLength = 2;

    quiver3(pos(1), pos(2), pos(3), R(1,1)*axisLength, R(2,1)*axisLength, R(3,1)*axisLength, 'color', 'r', 'LineWidth', 2);
    quiver3(pos(1), pos(2), pos(3), R(1,2)*axisLength, R(2,2)*axisLength, R(3,2)*axisLength, 'color', 'g', 'LineWidth', 2);
    quiver3(pos(1), pos(2), pos(3), R(1,3)*axisLength, R(2,3)*axisLength, R(3,3)*axisLength, 'color', 'b', 'LineWidth', 2);
end