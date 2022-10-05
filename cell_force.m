function f = cell_force(xvec,border,k1,k2,a0)
    kvec = zeros(1,length(xvec)-1);
    kvec(1,1:border-1) = k1;
    kvec(1,border:end) = k2;
    q = rolling_average_density(xvec);
    q1 = (1./q - a0);
    f = kvec.*q1;
end



