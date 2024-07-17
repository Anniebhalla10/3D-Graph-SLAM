function plotFunction(landmarksg, landmarks, landmarks_est, posesg, poses, poses_est)
    plotMesh();
    plotTrajectories(posesg, poses, poses_est);
    saveas(gcf,'Trajectories.png')
    
    plotMesh();
    plotMap(landmarksg, landmarks, landmarks_est);
    saveas(gcf,'Maps.png')
    
    plotMesh();
    plotTrajectories(posesg, poses, poses_est);
    plotMap(landmarksg, landmarks, landmarks_est);
    
    plot_all_landmarks(landmarksg);
    plot_all_poses(posesg);
    plot_errors();
end

function plotMesh()
    [xs, ys, zs] = sphere;
    R = 10;
    xs = R * xs;
    ys = R * ys;
    zs = R * zs;
    figure;
    mesh(xs, ys, zs, 'EdgeColor', [0.5, 0.5, 0.5]); hold on; axis equal;
end

function plotTrajectories(posesg, poses, poses_est)
    for i = 1:size(posesg,2)
        plotPose(posesg(:, i), 'r'); % Ground truth poses in red
    end
    for i = 1:size(poses,2)
        plotPose(poses(:, i), 'b'); % Initial guess poses in blue
    end
    for i = 1:size(poses_est,2)
        plotPose(poses_est(:, i), 'g'); % Estimated poses in green
    end
end

% Plot Pose
function hnd = plotPose(pose, color)
    hold on;
    
    pos = pose(1:3);
    R = Rx(pose(4))* Ry(pose(5))* Rz(pose(6)) ;
    
    axisLength = 2;
    
    hnd(1)=quiver3(pos(1), pos(2), pos(3), R(1,1)*axisLength, R(2,1)*axisLength, R(3,1)*axisLength, 'color', color, 'LineWidth', 2);
    hnd(2)=quiver3(pos(1), pos(2), pos(3), R(1,2)*axisLength, R(2,2)*axisLength, R(3,2)*axisLength, 'color', color, 'LineWidth', 2);
    hnd(3)=quiver3(pos(1), pos(2), pos(3), R(1,3)*axisLength, R(2,3)*axisLength, R(3,3)*axisLength, 'color', color, 'LineWidth', 2);
end

function plotMap(landmarksg, landmarks, landmarks_est)
    % Plot landmarks
    h1= scatter3(landmarksg(1,:), landmarksg(2,:), landmarksg(3,:), 'ro' ); % Ground truth landmarks in red
    h2= scatter3(landmarks(1,:), landmarks(2,:), landmarks(3,:), 'b', 'filled'); % Initial guess landmarks in blue
    h3 = scatter3(landmarks_est(1,:), landmarks_est(2,:), landmarks_est(3,:), 'g', 'filled'); % Estimated landmarks in green
    legend([h1, h2, h3], 'Ground Truth', 'Initial Guess', 'After LS');
end

function plot_all_landmarks(landmarksg)
    global all_landmarks;
    % Plot the mesh
    plotMesh();
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    scatter3(landmarksg(1,:), landmarksg(2,:), landmarksg(3,:), 'red','filled' ); % Ground truth landmarks in red
    
    % Plot all poses and landmarks
    for i = 1:length(all_landmarks)
        landmarks = all_landmarks{i};
        title('Optimization Evolution - Landmarks : iteration ',{i});
        
        if i>1
            delete(h1)
        end
        
        % Plot landmarks
        h1 = scatter3(landmarks(1,:), landmarks(2,:), landmarks(3,:), 'bo' );
        
        pause(0.5);
    end

end


function plot_all_poses(posesg)
    global all_poses;
    
    % Plot the mesh
    plotMesh();
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    for i = 1:size(posesg,2)
        plotPose(posesg(:, i), 'r'); % Ground truth poses in red
    end
    
    plots = [];
    
    % Plot all poses
    for i = 1:length(all_poses)
        poses = all_poses{i};
        title('Optimization Evolution - Poses : Iteration',{i});
        
        % Delete previous plots if they exist
        if ~isempty(plots)
            delete(plots(:));
        end
        
        % Plot poses
        num_poses = size(poses, 2);
        plots = gobjects(1, num_poses * 3); % 3 quivers per pose
        for j = 1:num_poses
            pose_plots = plotPose(poses(:, j), 'blue'); % Estimated poses in blue
            plots((j-1)*3+1:j*3) = pose_plots; % Store the plot handles
        end
        
        pause(0.5); % Adjust as needed
    end
end

function plot_errors()
    global chi_values;
    figure;
    plot(1:length(chi_values), chi_values, '-o', 'LineWidth', 2);
    
    % Formatting
    xlabel('Iteration');
    ylabel('\chi Value');
    title('Chi Values Evolution with Number of Iterations');
    grid on;
end

