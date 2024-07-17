function [poses, landmarks] = least_squares(landmarks, poses, transitions, observations, niterations, posesg)
    global all_poses all_landmarks chi_values;
    allp = numel(poses);

    % dim of all poses + landmarks 
    dim = allp + numel((landmarks));
    
    poses(:, 1) = posesg; 

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
            [e, Jr, Jl] = errorAndJacobianEuclideanLP(landmarks(:, idl), poses(:,pose_id), z);
                
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
            [e, Ji, Jj] = errorAndJacobianEuclideanPP(poses(:,pi_id), poses(:,pj_id), z);
            
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
function [e, Jr, Jl ] =errorAndJacobianEuclideanLP(XL, XR, z)
    t= XR(1:3);
    rx = Rx(XR(4));
    ry = Ry(XR(5));
    rz = Rz(XR(6));
    rx_p= Rx_prime(XR(4));
    ry_p = Ry_prime(XR(5));
    rz_p = Rz_prime(XR(6));

    z_hat= (rx*ry*rz)' *(XL - t);
    e = z_hat- z;
    Jr = zeros(3,6);
    Jr(1:3,1:3)= -(rx*ry*rz)';
    Jr(:,4)= (rx_p * ry* rz)' * (XL-t);
    Jr(:,5) =(rx * ry_p* rz)' * (XL-t);
    Jr(:,6)=(rx * ry* rz_p)' * (XL-t);
    Jl= (rx*ry*rz)';
end

% 2. Jacobian calculation for Pose Pose
function [e,Ji,Jj]= errorAndJacobianEuclideanPP(pi, pj, Z)
    % measurement function
    ti = pi(1:3);
    Ri = Rx(pi(4))*Ry(pi(5))*Rz(pi(6));
    tj = pj(1:3);
    Rj = Rx(pj(4))*Ry(pj(5))*Rz(pj(6));
    
    Ria = Rx_prime(pi(4))*Ry(pi(5))*Rz(pi(6));
    Rib = Rx(pi(4))*Ry_prime(pi(5))*Rz(pi(6));
    Ric = Rx(pi(4))*Ry(pi(5))*Rz_prime(pi(6));
    
    % tj as seen from pi for the translation part
    % for the rotation part,just the difference between the angle vectors
    Z_pred = [Ri'*(tj-ti); 
              pj(4:6)-pi(4:6)];
    % Error
    e = Z_pred -Z;

    I = eye(3,3);
    Ji = [-Ri' Ria'*(tj-ti) Rib'*(tj-ti) Ric'*(tj-ti);
          zeros(3,3), -I];
    Jj = [Ri' zeros(3,3);
          zeros(3,3) -I]; 
end

% 3. Update Poses and Landmarks
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
        np(:,i) = poses(:,i) + dxi;
    end
    % update landmarks
    nl = landmarks + reshape(updatel, 3, []);
end

% 4. Helper Scripts
function r_idx = poseMatrixIndex(pose_id)
    r_idx = 1+(pose_id-1)*6;
end

function l_idx = landmarkMatrixIndex(landmark_id, allp)
    l_idx = 1+ allp + (landmark_id-1)*3;
end

% 5.Calculate Rotation matrices derivatives
function rxp = Rx_prime(alpha)
       rxp = [0 0 0 ;
        0 -sin(alpha) -cos(alpha);
        0 cos(alpha) -sin(alpha)];
end

function ryp = Ry_prime(beta)
       ryp = [-sin(beta) 0 cos(beta) ;
        0 0 0;
        -cos(beta) 0 -sin(beta)];
end

function rzp = Rz_prime(gamma)
       rzp = [-sin(gamma) -cos(gamma) 0 ;
        cos(gamma) -sin(gamma) 0;
        0 0 0];
end
