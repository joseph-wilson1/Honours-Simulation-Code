function [xvec,kvec,avec,dead,N,border,L_index] = ...
    stochastic_birth(xvec,length_tol,b,kvec,avec,dead,multiple_cell,border,t)    
    % Function to implement birth mechanic.
    N = length(xvec)-1;

    [Lt,L_index] = L_t(xvec,length_tol); %Number of and vector of indices of 
                                          %cells that exceed minimum birth
                                          %length
    r = rand(Lt,1); %Vector of uniform random numbers.
    birth = nonzeros((r<b).*L_index'); %Cell indices that split.
    
    %Implement splitting of cells
    for k=1:length(birth)
        i = birth(k);
        %Dead cells cannot split
        if dead(i)
            continue
        end
        strcat("Birth at index = ",num2str(i), "at time", num2str(t), 'with length = ',num2str(li(xvec,birth(k))));
        %Split cell and create two copies.
        midpoint =  (xvec(i) + xvec(i+1))/2;
        %Introduce first copy at left x-value, and second copy at
        %midpoint x-value. Reassign x-values.
        xvec =  [xvec(1:i) midpoint xvec(i+1:N+1)]; %this is now N+2 long
        kvec = [kvec(1:i) kvec(i) kvec(i+1:N)]; %k_i and a_i are copied.
        avec = [avec(1:i) avec(i) avec(i+1:N)];
        dead = [dead(1:i) 0 dead(i+1:N)];
        N = N+1; %Increase number of cells count
        if multiple_cell %If cell birth occurs to left of border, shift
                         %border coord right, and reassign populations to
                           %correct space
            if i< border
                border=border+1;
            end
        end
        birth=birth+1; %Cell birth shifts cell index
        
        %If cell birth occurs to right of border, adjusting number of
        %cells reassigns populations to correct space.
    end
end