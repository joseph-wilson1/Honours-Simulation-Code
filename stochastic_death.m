function [kvec,avec,dead] = stochastic_death(N,d,kvec,avec,dead,k_death)

    %Stochastic Process - Death
    cells = 1:N;
    r = rand(N,1);
    death = nonzeros((r<d).*cells');
    
    %Implement death of cells
    for k=1:length(death)
        i = death(k);
        if dead(i)
            continue
        end
        avec(i) = 0; %Set cell length to zero.
        %kvec(i) = k_death*kvec(i); %Increase spring constant
        kvec(i) = k_death;
        dead(i) = 1;
        strcat('Cell death occurs for cell=',num2str(i));
    end

end