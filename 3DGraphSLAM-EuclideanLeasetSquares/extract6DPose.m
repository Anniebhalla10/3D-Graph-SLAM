function pose = extract6DPose(T)
  pos = T(1:3,4)';
  R = T(1:3, 1:3);
  
  % Handle gimbal lock 
  if abs(R(1,3)) > 0.95
    % Singularity case (z-axis)
    rot_x = 0;
    rot_y = pi/2;
    rot_z = 0;
    % pos(1:3) = zeros(1,3);
  else
    % No singularity
    rot_x = atan2(-R(2,3), R(3,3));
    rot_y = asin(R(1,3));
    rot_z = atan2(-R(1,2), R(1,1));
  end

  rot_x = mod(rot_x + pi, 2*pi) - pi;
  rot_y = mod(rot_y + pi, 2*pi) - pi;
  rot_z = mod(rot_z + pi, 2*pi) - pi;
  
  % Create a 6D vector with position and Euler angles
  pose = [pos rot_x rot_y rot_z]';
end