function mvec = midp_cell(xvec)
    N=length(xvec);
    mvec=zeros(1,N-1);
    for i=1:N-1
        mvec(i) = (xvec(i)+xvec(i+1))/2;
    end
end