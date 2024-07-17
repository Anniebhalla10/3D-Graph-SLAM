% v = [x,y,z, alpha,beta,gamma]
% A = 4*4 homogeneous transformation matrix
function A = v2t(v)
    x= v(1);
    y= v(2);
    z= v(3);
    alpha= v(4);
    beta= v(5);
    gamma= v(6);

    R_x= [1,0,0;
        0, cos(alpha), -sin(alpha);
        0, sin(alpha), cos(alpha)];

    R_y = [cos(beta),0,sin(beta);
        0, 1, 0;
        -sin(beta), 0, cos(beta)];

     R_z = [cos(gamma),-sin(gamma),0;
        sin(gamma), cos(gamma), 0;
        0, 0, 1];

     R = R_x * R_y* R_z;

     A = [R, [x;y;z]
         0,0,0,1];
end