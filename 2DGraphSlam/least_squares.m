function [poses, landmarks] = least_square(landmarks, poses, transitions, observations, id_to_landmark, niterations, posesg)
    global all_poses all_landmarks
    % dim of all poses (imagine a column vector with all poses stacked) 
    allp = numel((poses));
    % dim of all poses + landmarks 
    dim = allp + numel((landmarks));
    % NOW WE FIX 3 POINTS SINCE THE PROBLEM AS IT IS HAS 3DOF (from 2d rigid transform)
    poses(:, 1:3) = posesg; 
    for iteration = 1:niterations
        H = zeros(dim, dim);  b = zeros(dim, 1); %accumulators for H and b
        chi_tot = 0; %cumulative chi2
        % iterations over all observations (landmark pose optimization)
        for i = 1:length(observations)
            % current observation
            obs = observations(i);
            % pose from where we observe
            pose_id = obs.pose_id;
            % iterations over each observed landmark
            for j = 1:length(obs.observation)
                z = obs.observation(j);
                % id of the observed landmark
                id = z.id;
                % measurement for the observed landmark
                d = z.measurement;
                % we retrieve the landmark vector position from its id
                idl = id_to_landmark(id);
                % Gauss Newton algorithm
                [e, Jr, Jl] = errorAndJacobianLP(landmarks(:, idl), poses(:, pose_id), d);
                chi = e' * e;
                %  simplest robust kernel
                k_treshold = 3;
                if chi > k_treshold
                    e = e*sqrt(k_treshold/chi);
                    chi = k_treshold;
                end
                chi_tot = chi_tot + chi;
                %  now we need to compute H and b
                % index for contribution of the current pose in the H and b matrices
                rhidx = 1 + (pose_id-1) * 3;
                % index for contribution of the current landmark in the H and b matrices
                lhidx = 1 + allp + (idl-1) * 2;
                % product of the jacobians that appear in the H matrix
                Jrl = Jr' * Jl;
                % contribution to the H matrix
                H(rhidx:rhidx+2, rhidx:rhidx+2) = H(rhidx:rhidx+2, rhidx:rhidx+2) + Jr' * Jr;
                H(lhidx:lhidx+1, lhidx:lhidx+1) = H(lhidx:lhidx+1, lhidx:lhidx+1) + Jl' * Jl;
                H(rhidx:rhidx+2, lhidx:lhidx+1) = H(rhidx:rhidx+2, lhidx:lhidx+1) + Jrl;
                H(lhidx:lhidx+1, rhidx:rhidx+2) = H(lhidx:lhidx+1, rhidx:rhidx+2) + Jrl';
                % contribution to the b matrix
                b(rhidx:rhidx+2) = b(rhidx:rhidx+2) + Jr' * e;
                b(lhidx:lhidx+1) = b(lhidx:lhidx+1) + Jl' * e;              
            end
        end 
        % iterations over all transitions (pose pose optimization)
        for i = 1:length(transitions)
            % current transition
            trans = transitions(i);
            % pose from where we translate
            pi_id = trans.id_from;
            % pose to where we translate
            pj_id = trans.id_to;
            % odometry measurement
            z = trans.v;
            % Gauss Newton algorithm (with a simple robust kernel)
            [e, Ji, Jj] = errorAndJacobianPP(poses(:, pi_id), poses(:, pj_id), z);
            k_treshold = 1.5;
            if chi > k_treshold
                e = e*sqrt(k_treshold/chi);
                chi = k_treshold;
            end
            chi_tot = chi_tot + chi;
            % index for contribution of the pose i in the H and b matrices
            iidx = 1 + (pi_id - 1) * 3;
            % index for contribution of the pose j in the H and b matrices
            jidx = 1 + (pj_id - 1) * 3;
            % product of the jacobians that appear in the H matrix
            jij = Ji' * Jj;
            % contribution to the H matrix
            H(iidx:iidx+2, iidx:iidx+2) = H(iidx:iidx+2, iidx:iidx+2) + Ji' * Ji;
            H(jidx:jidx+2, jidx:jidx+2) = H(jidx:jidx+2, jidx:jidx+2) + Jj' * Jj;
            H(iidx:iidx+2, jidx:jidx+2) = H(iidx:iidx+2, jidx:jidx+2) + jij;
            H(jidx:jidx+2, iidx:iidx+2) = H(jidx:jidx+2, iidx:iidx+2) + jij';
            % contribution to the b matrix
            b(iidx:iidx+2) = b(iidx:iidx+2) + Ji' * e;
            b(jidx:jidx+2) = b(jidx:jidx+2) + Jj' * e;   
        end
        % initialize dx vector (the fixed poses are not optimized)
        dx = zeros(dim, 1);
        H = H(4:end, 4:end);
        b = b(4:end);
        tdx = -H\b;
        dx(4:end) = tdx;
        % update poses and landmarks with the computed dx
        [poses, landmarks] = update(poses, landmarks, dx, allp);
        fprintf("chi after iteration %d: %f\n", iteration, chi_tot)
        all_poses{iteration+1} = poses;
        all_landmarks{iteration+1} = landmarks;
    end        
end

function [e, Jr, Jl] = errorAndJacobianLP(land, pos, z)
    %  standard jacobian for 2D obs
    Xr = v2t(pos);
    R = Xr(1:2, 1:2);
    t = Xr(1:2, 3);
    zh = R' * (land - t);
    e = zh - z;   
    Jr = zeros(2,3);
    Jr(1:2,1:2) = -R';
    Jr(1:2,3) = R' * [0 1; -1 0] * (land - t);
    Jl = R';
end

function [e, Ji, Jj] = errorAndJacobianPP(pi, pj, z)
    % here we apply the compounding operators for odometry error and Jacobian
    xij = CF(pi,pj);
    e = xij - z;
    % always remember to wrap the angle
    e(3)=wrapToPi(e(3));
    Jij = JF(pi,pj);
    [Ji, Jj] = deal(Jij(:,1:3), Jij(:,4:6));
end

function [np, nl] = update(poses, landmarks, dx, allp)
    % split dx in pose and landmark updates
    updatep = dx(1:allp);
    updatel = dx(allp+1:end);
    % update poses
    upvp = reshape(updatep, 3, []);
    np = poses;
    for i = 1:size(upvp,2)
        dxi = upvp(:, i);
        % compounding operator to add the dx
        np(:, i) = CP(dxi,poses(:, i));
    end
    % update landmarks
    nl = landmarks + reshape(updatel, 2, []);
end

% computes the homogeneous transform matrix A of the pose vector v
% A:[ R t ] 3x3 homogeneous transformation matrix, r translation vector
% v: [x,y,theta]  2D pose vector
function A = v2t(v)
  	c = cos(v(3));
  	s = sin(v(3));
	A = [c, -s, v(1);
         s, c, v(2);
         0, 0, 1];
end

%% functions from the paper: R. Smith et al., Estimating Uncertain Spatial Relationships in Robotics
% compounding operator
function xik = CP(xij,xjk)
    xik = [xjk(1)*cos(xij(3))-xjk(2)*sin(xij(3)) + xij(1); ...
        xjk(1)*sin(xij(3))+xjk(2)*cos(xij(3)) + xij(2); ...
        xij(3) + xjk(3)];
end

% jacobian of compounding operator
function J = JP(xij,xjk)
    xik = CP(xij,xjk);
    J = [ 1 0 -(xik(2)-xij(2)) cos(xij(3)) -sin(xij(3)) 0; ...
        0 1 (xik(1)-xij(1)) sin(xij(3)) cos(xij(3)) 0; ...
        0 0 1 0 0 1];
end

% inverse compounding operator
function xik = CI(xij)
    xik = [-xij(1)*cos(xij(3))-xij(2)*sin(xij(3)); ...
        xij(1)*sin(xij(3))-xij(2)*cos(xij(3)); ...
        -xij(3)];
end

% jacobian of inverse compounding operator
function J = JI(xij)
    xji = CI(xij);
    J = [-cos(xij(3)) -sin(xij(3)) xji(2); ...
        sin(xij(3)) -cos(xij(3)) -xji(1); ...
        0 0 -1];
end

% composition for xwi, xwj -> xij
function xij = CF(xwi,xwj)
    xij = CP(CI(xwi),xwj);
end

% jacobian for the composition of compounding operators above
function J = JF(xwi,xwj)
    xiw = CI(xwi);
    jpij = JP(xiw,xwj);
    jm = JI(xwi);
    J = [jpij(:,1:3)*jm jpij(:,4:6)];
end
