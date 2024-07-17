function X_inv = inverse(X)
    R = X(1:3,1:3);
    t = X(1:3,4);

    R_inv = R';
    t_inv= -R_inv * t;

    X_inv = [R_inv, t_inv;
        0,0,0,1];
end