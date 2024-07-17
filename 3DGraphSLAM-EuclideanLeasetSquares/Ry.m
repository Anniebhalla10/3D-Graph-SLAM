function R_Y = Ry(beta)
    R_Y = [cos(beta) 0 sin(beta) ;
        0 1 0;
        -sin(beta) 0 cos(beta)];
end