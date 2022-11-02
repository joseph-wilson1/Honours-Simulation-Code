function [xvec_new,kvec,avec,N,border] = stochastic_death2(N,d,x,k,a,border)
    % Function to implement instant death mechanic in each time step.

    x_old = x;
    k_old = k;
    a_old = a;
    tol = 1e-4;

    %Stochastic Process - Death
    cells = 1:N;
    r = rand(N,1);
    death = nonzeros((r<d).*cells');
    
    %Implement death of cells
    for k=1:length(death)
        i = death(k);
        if abs(i - 1) < tol
            x_old = [x_old(1) x_old(3:end)];
            k_old = k_old(2:end);
            a_old = a_old(2:end);
            N=N-1;
            border = border-1;
        elseif abs(i - N) < tol
            x_old = [x_old(1:N-1) x_old(end)];
            k_old = k_old(1:end-1);
            a_old = a_old(1:end-1);
            N = N-1;
        else
            li = x_old(i+1) - x_old(i);
            xi_new = x(i) + li/2;
            x_old = [x_old(1:i-1) xi_new x_old(i+2:end)];
            k_old = [k_old(1:i-1) k_old(i+1:end)];
            a_old = [a_old(1:i-1) a_old(i+1:end)];
            N = N-1;
            if i < border
                border=border-1;
            end
        end
        strcat('Cell death occurs for cell=',num2str(i));
        death=death-1;
    end

    xvec_new = x_old;
    kvec = k_old;
    avec = a_old;
end