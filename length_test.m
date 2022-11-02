function [xvec_new,kvec,avec,deadvec,N,border] = length_test(x,k,a,dead,l_tol,border)
%Function that discerns wether or not a cell is small enough to be removed,
%and then removes it from the simulation, along with its attributes. For
%gradual death option.

N = length(x)-1; %Number of cells.
x_old = x;
k_old = k;
a_old = a;
dead_old = dead;
if li(x_old,1) < l_tol
    x_old = [x_old(1) x_old(3:end)];
    k_old = k_old(2:end);
    a_old = a_old(2:end);
    dead_old = dead_old(2:end);
    N=N-1;
    border = border-1;
elseif li(x_old,N) < l_tol
    x_old = [x_old(1:N-1) x_old(end)];
    k_old = k_old(1:end-1);
    a_old = a_old(1:end-1);
    dead_old = dead_old(1:end-1);
    N=N-1;
end
i = 2;
while i<N
    if li(x_old,i) < l_tol
        xlen = x_old(i+1) - x_old(i);
        xi_new = x(i) + xlen/2;
        x_old = [x_old(1:i-1) xi_new x_old(i+2:end)];
        k_old = [k_old(1:i-1) k_old(i+1:end)];
        a_old = [a_old(1:i-1) a_old(i+1:end)];
        dead_old = [dead_old(1:i-1) dead_old(i+1:end)];
        N = N-1;
        if i < border
            border=border-1;
        end
    else
        i=i+1;
    end
end
xvec_new = x_old;
kvec = k_old;
avec = a_old;
deadvec = dead_old;
end