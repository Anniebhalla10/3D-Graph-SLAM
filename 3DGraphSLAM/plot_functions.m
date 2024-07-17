function plot_functions(landmarks, landmarksg, landmarks_est, poses, posesg, poses_est)
    % plot of the trajectories
    plot_trajectories(poses(1:3,4,:), poses_est(1:3,4,:), posesg(1:3,4,:))
    saveas(gcf,'Trajectories.png')
   
    % plot of the maps (landmarks)
    plot_maps(landmarks, landmarks_est, landmarksg, false)
    saveas(gcf,'Landmarks.png')
    
    % both the trajectories and the maps are plotted in the same figure
    plot_trajectories(poses(1:3,4,:), poses_est(1:3,4,:), posesg(1:3,4,:))
    plot_maps(landmarks, landmarks_est, landmarksg, true)
    saveas(gcf,'TrajectoriesAndLandmarks.png')  
end

function plot_maps(landmarks, landmarksest, landmarksg, holdon)
    if holdon
        hold on;
        title("Full result");
    else
        figure();
        hold on;
        axis equal;
        title("Maps comparison");
    end

    % Plotting landmarks
     plot3(landmarks(1, :), landmarks(2, :), landmarks(3, :), 'o', 'MarkerEdgeColor', 'blue', 'MarkerSize', 10, 'LineWidth', 2);
     plot3(landmarksg(1, :), landmarksg(2, :), landmarksg(3, :), 'o', 'MarkerEdgeColor', 'red', 'MarkerSize', 10, 'LineWidth', 2);
     plot3(landmarksest(1, :), landmarksest(2, :), landmarksest(3, :), 'o', 'MarkerEdgeColor', 'green', 'MarkerSize', 10, 'LineWidth', 4);
   
    % Legend
    legend("Initial Guess", "Ground Truth","After LS");

    % View
    view(3);

    % Additional Formatting (optional)
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');

    % Hold off if new figure was created
    if ~holdon
        hold off;
    end
end


  function plot_trajectories(poses, posesest, posesg)
    num_points = size(poses, 3);

    figure();
    hold on;
    plot3(reshape(poses(1, 1, :), 1, num_points), ...
          reshape(poses(2, 1, :), 1, num_points), ...
          reshape(poses(3, 1, :), 1, num_points), 'b-', 'LineWidth', 2);
      
    plot3(reshape(posesg(1, 1, :), 1, num_points), ...
          reshape(posesg(2, 1, :), 1, num_points), ...
          reshape(posesg(3, 1, :), 1, num_points), 'r-', 'LineWidth', 2);

    plot3(reshape(posesest(1, 1, :), 1, num_points), ...
          reshape(posesest(2, 1, :), 1, num_points), ...
          reshape(posesest(3, 1, :), 1, num_points), 'g-', 'LineWidth', 2);
      
    % Formatting
    axis equal;
    title('Trajectories Comparison');
    legend('Initial Guess', 'Ground Truth','After LS');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    view(3); 

    hold off;
  end
