% Returns a 3Dtensor for poses where each ith index in the third dimension
% has a homogeneous transofrmation of the robot poses
% Robot poses are coordinates of the robots as seen from the world

function poses_ht = create_HT_poses(poses)
    num_poses= size(poses,2);
    poses_ht = zeros(4,4,num_poses);    
        for i = 1:num_poses
            poses_ht(:,:,i) = v2t(poses(:,i));
        end
end