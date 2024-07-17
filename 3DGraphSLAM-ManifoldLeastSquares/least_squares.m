function [poses, landmarks] = least_squares(landmarks, poses, transitions, observations, niterations, posesg)
    global all_landmarks all_poses chi_values
    allp = size(poses,3)*6;

    % dim of all poses + landmarks 
    dim = allp + numel((landmarks));
    
    poses(:, :, 1) = posesg; 

    for iteration= 1:niterations
        H = zeros(dim, dim);  b = zeros(dim, 1); %accumulators for H and b
        chi_tot = 0; %cumulative chi2

        % 1. Landmark Pose Optimization
         % iterations over all observations (landmark pose optimization)
        for i = 1:length(observations)
            % current observation
            obs = observations(i);
            % pose from where we observe
            pose_id = obs.p_idx;
            idl= obs.l_idx;
            z = obs.meas;
            % Gauss Newton algorithm
            [e, Jr, Jl] = errorAndJacobianManifoldLP(landmarks(:, idl), poses(:, :,pose_id), z);
                
            chi = e' * e;
            k_treshold = 1.2;
            if chi > k_treshold
               e = e * sqrt(k_treshold/chi);
               chi = k_treshold;
            end
            chi_tot = chi_tot + chi;
                
            % Computing H and b
            % index for contribution of the current pose in the H and b matrices
            rhidx = poseMatrixIndex(pose_id);
            %index for contribution of the current landmark in the H and b matrices
            lhidx = landmarkMatrixIndex(idl, allp);
                 
            % product of the jacobians that appear in the H matrix
            Jrl = Jr' * Jl;
            % contribution to the H matrix
            H(rhidx:rhidx+5, rhidx:rhidx+5) = H(rhidx:rhidx+5, rhidx:rhidx+5) + Jr' * Jr;
            H(lhidx:lhidx+2, lhidx:lhidx+2) = H(lhidx:lhidx+2, lhidx:lhidx+2) + Jl' * Jl;
            H(rhidx:rhidx+5, lhidx:lhidx+2) = H(rhidx:rhidx+5, lhidx:lhidx+2) + Jrl;
            H(lhidx:lhidx+2, rhidx:rhidx+5) = H(lhidx:lhidx+2, rhidx:rhidx+5) + Jrl';
            %contribution to the b matrix
            b(rhidx:rhidx+5) = b(rhidx:rhidx+5) + Jr' * e;
            b(lhidx:lhidx+2) = b(lhidx:lhidx+2) + Jl' * e;              
        end 

        % 2. Pose-Pose Optimization
        for i = 1:length(transitions)
            % current transition in homogeneos coordinates
            trans = transitions(i);
          
            % pose from where we translate
            pi_id = trans.from;
            % pose to where we translate
            pj_id = trans.to;
            % odometry measurement
            z = trans.meas;
       
            % Gauss Newton algorithm (with a simple robust kernel)
            [e, Ji, Jj] = errorAndJacobianManifoldPP(poses(:, :,pi_id), poses(:, :,pj_id), z);
            
            chi= e' *e;
            k_treshold = 1.5;
            if chi > k_treshold
                e = e*sqrt(k_treshold/chi);
                chi = k_treshold;
            end
            chi_tot = chi_tot + chi;

            % index for contribution of the pose i in the H and b matrices
            iidx = poseMatrixIndex(pi_id);
            % index for contribution of the pose j in the H and b matrices
            jidx = poseMatrixIndex(pj_id);
            % product of the jacobians that appear in the H matrix
            Jij = Ji' * Jj;
            % contribution to the H matrix
            H(iidx:iidx+5, iidx:iidx+5) = H(iidx:iidx+5, iidx:iidx+5) + Ji' * Ji;
            H(jidx:jidx+5, jidx:jidx+5) = H(jidx:jidx+5, jidx:jidx+5) + Jj' * Jj;
            H(iidx:iidx+5, jidx:jidx+5) = H(iidx:iidx+5, jidx:jidx+5) + Jij;
            H(jidx:jidx+5, iidx:iidx+5) = H(jidx:jidx+5, iidx:iidx+5) + Jij';
            % contribution to the b matrix
            b(iidx:iidx+5) = b(iidx:iidx+5) + Ji' * e;            
            b(jidx:jidx+5) = b(jidx:jidx+5) + Jj' * e;   
        end

        % initialize dx vector (the fixed poses are not optimized)
        dx = zeros(dim, 1);
        H = H(7:end, 7:end);
        b = b(7:end);
        tdx = -(H+1e-6*eye(size(H)))\b;
        dx(7:end) = tdx;
        % update poses and landmarks with the computed dx
        [poses, landmarks] = update(poses, landmarks, dx, allp);
        fprintf("chi after iteration %d: %f\n", iteration, chi_tot)
        chi_values(iteration+1)= chi_tot;
        all_poses{iteration+1} = poses;
        all_landmarks{iteration+1} = landmarks;
    end
end


% 1. Jacobian calculation for Landmark Pose
function [e, Jr, Jl ] =errorAndJacobianManifoldLP(XL, XR, z)
    R = XR(1:3,1:3);
    t= XR(1:3,4);
    z_hat= R' *(XL - t);
    e = z_hat- z;
    Jr = zeros(3,6);
    Jr(1:3,1:3)= eye(3);
    Jr(1:3, 4:6)= -skew(XL);
    Jr = -R'*Jr;
    Jl= R';
end

% 2. Jacobian calculation for Pose Pose
function [e,Ji,Jj]= errorAndJacobianManifoldPP(pi, pj, Z)
    Rdx= [0,0,0;
        0,0,-1;
        0,1,0];

    Rdy= [0,0,1;
        0,0,0;
        -1,0,0];

    Rdz= [0,-1,0;
        1,0,0;
        0,0,0];

    Xi= pi;
    Xi_inv = inverse(Xi);
    Ri_transpose = Xi_inv(1:3,1:3); 

    Rj = pj(1:3,1:3);
    tj = pj(1:3,4);

    rxp = Ri_transpose* Rdx* Rj;
    rxp = reshape(rxp,[9 1]);

    ryp = Ri_transpose*Rdy * Rj;
    ryp = reshape(ryp,[9 1]);

    rzp = Ri_transpose * Rdz * Rj;
    rzp = reshape(rzp,[9 1]);

    Jj = [zeros(9,3) , rxp, ryp, rzp;
          Ri_transpose, -Ri_transpose*skew(tj)];
    Ji = -Jj;

    g= Xi_inv*pj;
  
    e =  flatten(g) - flatten(Z);
end

% 3. Update Poses
function [np, nl] = update(poses, landmarks, dx, allp)
    % split dx in pose and landmark updates
    updatep = dx(1:allp);
    updatel = dx(allp+1:end);
    % update poses
    upvp = reshape(updatep, 6, []);
    np = poses;
    for i = 1:size(upvp,2)
        dxi = upvp(:, i);
        % box plus operator to add the dx to the initial guess
        np(:,:, i) = boxPlus_pose(poses(:,:,i),dxi);
    end
    % update landmarks
    nl = landmarks + reshape(updatel, 3, []);
end

% 4. Box Plus
% XR = robot poses
% XL = landmark poses
function XR = boxPlus_pose(Xr, dx_r)
    % Adding dx to robot pose
        XR = v2t(dx_r)* Xr;
end


% Helper Functions

function X =  flatten(transition)
    X = transition(1:size(transition,1)-1,:);
    X = reshape(X, [12,1]);
end

function r_idx = poseMatrixIndex(pose_id)
    r_idx = 1+(pose_id-1)*6;
end

function l_idx = landmarkMatrixIndex(landmark_id, allp)
    l_idx = 1+ allp + (landmark_id-1)*3;
end
