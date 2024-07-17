function R_X = Rx(alpha)
    R_X = [1 0 0 ;
        0 cos(alpha) -sin(alpha);
        0 sin(alpha) cos(alpha)];
end