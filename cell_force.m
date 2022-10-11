function f = cell_force(xvec,border,k1,k2,a0)
    kvec = zeros(1,length(xvec)-1);
    kvec(1,1:border-1) = k1;
    kvec(1,border:end) = k2;
    q = density_vector(xvec);
    q1 = (a0 - 1./q);
    f = kvec.*q1;
end



