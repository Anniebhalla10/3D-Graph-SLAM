function R_Z = Rz(gamma)
    R_Z = [cos(gamma) -sin(gamma) 0 ;
        sin(gamma) cos(gamma) 0;
        0 0 1];
end