function [final_time1,final_time2,density,border_pos] = cell_invasion(k1,k2,L,N,n,wall_n,b,d,dt)

    %n is simulation steps, T = n*dt
    n_initial = wall_n; %Steps before wall is removed.
    eta = 1; %drag coefficient

    %d is death Probability 
    %b is birth Probability

    %Multiple cell flag
    multiple_cell = 1;

    x0 = 0; %Left boundary of domain
    k_hom = 0.5; %Spring coefficient (if homogenous)
    a0=L/N; %Resting cell length (homogenous)

    density = zeros(multiple_cell+1,n+1);%Vector of total density of each cell

    xvec = linspace(x0,L,N+1); %Set up spatial coordinates

    %nitial inhomogeneity of spatial positions of cells borders
    %xvec(2) = xvec(2) + 1/2 * L/N;
    %xvec(19) = xvec(19) + 1/2 * L/N;
    %xvec(22) = xvec(20);
    %xvec(10) = xvec(11) - 0.1*L/N;

    dead = zeros(1,N); %Vector that stores wether a cell is "dead" and reducing.
    death_count_1 = zeros(n+1,1); %Amount of dead cells of type 1 in time step j
    death_count_2 = zeros(n+1,1); %Amount of dead cells of type 2 in time step j
    kvec = zeros(1,N); %Spring coefficient vector for heterogenous population.
    a = zeros(1,N); %Resting cell length vector for heterogenous population.
    %Initial vector of resting cell lengths
    a(1,:) = a0;

    if multiple_cell == 1
        border = round(N/2); %spatial coord of border between cell types, 
                             %initially at midpoint of cell populations.
        %Set up multiple cell populations
        kvec(1,1:border) = k1; 
        kvec(1,border+1:N) = k2;
        N1 = border;
        xvec1 = xvec(1:N1+1);
        kvec1 = kvec(1:N1);
        a1 = a(1:N1);
        dead1 = dead(1:N1);
        N2 = border;
        xvec2 = xvec(N1+1:end);
        kvec2 = kvec(N2+1:end);
        a2 = a(N2+1:end);
        dead2 = dead(N2+1:end);
        %Initial density is equal.
        density(1,1) = border;
        density(2,1) = N-border;
    else
        %Set up single cell population
        kvec(1,:) = k_hom;
        %Initial density
        density(1) = N;
    end

    %Birth length threshold
    l_birth = 1.1*a0;
    
    border_pos(1) = xvec(border);

    %What is this??
    kvec(10) = k_hom*1.1;

    length_tol = a0/10; %Length under which we remove cell.
    k_death = 4; %Factor increase at death of spring constant.

    T = n*dt;

    'Start'
    break_val = 0; %Flag to debug with.
    while break_val == 0    
        for j=1:n %simulation steps (time)
            if j <= n_initial
                %Run simulation for first population!!
                if li(xvec1,1) < length_tol && N1 > 1
                    strcat('Cell =',num2str(1),' removed with length = ',num2str(li(xvec1,1)));
                    xvec1 = [xvec1(1) xvec1(3:N1+1)];
                    kvec1 = kvec1(2:N1);
                    a1 = a1(2:N1);
                    dead1 = dead1(2:N1);
                    N1  = N1-1; %Reduce number of cells.
                elseif li(xvec1,1) < length_tol
                    N1 = 0; % Finishing rules (to add...)
                end
                i = 2; %cell-border index step (Loop from second cell to second last cell)
                while i<=N1 % (While loop,as N can change within loop)
                    %Equation of motion
                    xvec_new = xvec1(i) + hookfun(eta,fi(li(xvec1,i),kvec1(i), ...
                        a1(i)),fi(li(xvec1,i-1),kvec1(i-1),a1(i-1)))*dt;
                    %Catch overshooting
                    if xvec_new < xvec1(i-1)
                        xvec_new = xvec1(i-1)*1.001;
                    elseif xvec_new > xvec1(i+1)
                        xvec_new = xvec1(i+1) * 0.999;
                    end
                    xvec1(i) = xvec_new; %Reassign new position.

                    %If length of cell becomes too small, remove cell
                    if li(xvec1,i) < length_tol
                        death = 1;
                        strcat('Cell =',num2str(i),' removed with length = ',num2str(li(xvec1,i)));
                        xvec1 = [xvec1(1:i-1) xvec1(i+1:N1+1)];
                        kvec1 = [kvec1(1:i-1) kvec1(i+1:N1)];
                        a1 = [a1(1:i-1) a1(i+1:N1)];
                        dead1 = [dead1(1:i-1) dead1(i+1:N1)];
                        N1  = N1-1; %Reduce number of cells.
                    end
                    i = i+1; %Step cell-border index.
                end

                %Stochastic process - Birth
                [Lt,L_index] = L_t(xvec1,l_birth); %Number of and vector of indices of 
                                                  %cells that exceed minimum birth
                                                  %length
                r = rand(Lt,1); %Vector of uniform random numbers.
                birth = nonzeros((r<b).*L_index'); %Cell indices that split.

                %Implement splitting of cells
                for k=1:length(birth)
                    i = birth(k);
                    %Dead cells cannot split
                    if dead1(i)
                        continue
                    end
                    strcat("Birth at index = ",num2str(i), 'in pop 1 at time', num2str(j*dt));
                    %Split cell and create two copies.
                    midpoint =  (xvec1(i) + xvec1(i+1))/2;
                    %Introduce first copy at left x-value, and second copy at
                    %midpoint x-value. Reassign x-values.
                    xvec1 =  [xvec1(1:i) midpoint xvec1(i+1:N1+1)]; %this is now N+2 long
                    kvec1 = [kvec1(1:i) kvec1(i) kvec1(i+1:N1)]; %k_i and a_i are copied.
                    a1 = [a1(1:i) a1(i) a1(i+1:N1)];
                    dead1 = [dead1(1:i) 0 dead1(i+1:N1)];
                    N1 = N1+1; %Increase number of cells count
                end

                %Stochastic Process - Death
                cells = 1:N1;
                r = rand(N1,1);
                death = nonzeros((r<d).*cells');

                %Implement death of cells
                for k=1:length(death)
                    i = death(k);
                    a1(i) = 0; %Set cell length to zero.
                    kvec1(i) = k_death*kvec1(i); %Increase spring constant
                    dead1(i) = 1;
                    strcat('Cell death occurs for cell=',num2str(i), 'in pop 1 at time=', num2str(j*dt));
                end


                %Run simulation for SECOND population!!
                if li(xvec2,1) < length_tol && N2 > 1
                    strcat('Cell =',num2str(1),' removed with length = ',num2str(li(xvec2,1)));
                    xvec2 = [xvec2(1) xvec2(3:N2+1)];
                    kvec2 = kvec2(2:N2);
                    a2 = a2(2:N2);
                    dead2 = dead2(2:N2);
                    N2  = N2-1; %Reduce number of cells.
                elseif li(xvec2,1) < length_tol
                    N2 = 0; % Finishing rules (to add...)
                end
                i = 2; %cell-border index step (Loop from second cell to second last cell)
                while i<=N2 % (While loop,as N can change within loop)
                    %Equation of motion
                    xvec_new = xvec2(i) + hookfun(eta,fi(li(xvec2,i),kvec2(i), ...
                        a2(i)),fi(li(xvec2,i-1),kvec2(i-1),a2(i-1)))*dt;
                    %Catch overshooting
                    if xvec_new < xvec2(i-1)
                        xvec_new = xvec2(i-1)*1.001;
                    elseif xvec_new > xvec2(i+1)
                        xvec_new = xvec2(i+1) * 0.999;
                    end
                    xvec2(i) = xvec_new; %Reassign new position.

                    %If length of cell becomes too small, remove cell
                    if li(xvec2,i) < length_tol
                        death = 1;
                        strcat('Cell =',num2str(i),' removed with length = ',num2str(li(xvec2,i)));
                        xvec2 = [xvec2(1:i-1) xvec2(i+1:N2+1)];
                        kvec2 = [kvec2(1:i-1) kvec2(i+1:N2)];
                        a2 = [a2(1:i-1) a2(i+1:N2)];
                        dead2 = [dead2(1:i-1) dead2(i+1:N2)];
                        N2  = N2-1; %Reduce number of cells.
                    end
                    i = i+1; %Step cell-border index.
                end

                %Stochastic process - Birth
                [Lt,L_index] = L_t(xvec2,l_birth); %Number of and vector of indices of 
                                                  %cells that exceed minimum birth
                                                  %length
                r = rand(Lt,1); %Vector of uniform random numbers.
                birth = nonzeros((r<b).*L_index'); %Cell indices that split.

                %Implement splitting of cells
                for k=1:length(birth)
                    i = birth(k);
                    %Dead cells cannot split
                    if dead2(i)
                        continue
                    end
                    strcat("Birth at index = ",num2str(i), 'in pop 2 at time', num2str(j*dt)); 
                    %Split cell and create two copies.
                    midpoint =  (xvec2(i) + xvec2(i+1))/2;
                    %Introduce first copy at left x-value, and second copy at
                    %midpoint x-value. Reassign x-values.
                    xvec2 =  [xvec2(1:i) midpoint xvec2(i+1:N2+1)]; %this is now N+2 long
                    kvec2 = [kvec2(1:i) kvec2(i) kvec2(i+1:N2)]; %k_i and a_i are copied.
                    a2 = [a2(1:i) a2(i) a2(i+1:N2)];
                    dead2 = [dead2(1:i) 0 dead2(i+1:N2)];
                    N2 = N2+1; %Increase number of cells count
                end

                %Stochastic Process - Death
                cells = 1:N2;
                r = rand(N2,1);
                death = nonzeros((r<d).*cells');

                %Implement death of cells
                for k=1:length(death)
                    i = death(k);
                    a2(i) = 0; %Set cell length to zero.
                    kvec2(i) = k_death*kvec2(i); %Increase spring constant
                    dead2(i) = 1;
                    strcat('Cell death occurs for cell=',num2str(i), 'in pop 2 at time=', num2str(j*dt));
                end

                %Account of density after stochastic processes
                if multiple_cell
                    density(1,j+1) = N1;
                    density(2,j+1) = N2;
                else
                    density(j+1) = N;
                end

                %Account of dead cells
                if multiple_cell
                    death_count_1(j) = sum(dead1);
                    death_count_2(j) = sum(dead2);
                else
                    death_count_1(j) = sum(dead(:));
                end

                %Account of master variables after wall is removed
                if j == n_initial
                    N = N1 + N2;
                    border = N1;
                    xvec = [xvec1 xvec2(2:end)];
                    kvec = [kvec1 kvec2];
                    a = [a1 a2];
                    dead = [dead1 dead2];
                end
            else
                i = 2; %cell-border index step (Loop from second cell to second last cell)
                if li(xvec,1) < length_tol
                    strcat('Cell =',num2str(1),' removed with length = ',num2str(li(xvec,1)));
                    xvec = [xvec(1) xvec(3:N+1)];
                    kvec = kvec(2:N);
                    a = a(2:N);
                    dead = dead(2:N);
                    N  = N-1; %Reduce number of cells.
                    border = border-1;
                end
                while i<=N % (While loop,as N can change within loop)
                    %Equation of motion
                    xvec_new = xvec(i) + hookfun(eta,fi(li(xvec,i),kvec(i), ...
                        a(i)),fi(li(xvec,i-1),kvec(i-1),a(i-1)))*dt;
                    %Catch overshooting
                    if xvec_new < xvec(i-1)
                        xvec_new = xvec(i-1)*1.001;
                    elseif xvec_new > xvec(i+1)
                        xvec_new = xvec(i+1) * 0.999;
                    end
                    xvec(i) = xvec_new; %Reassign new position.

                    %If length of cell becomes too small, remove cell
                    if li(xvec,i) < length_tol
                        strcat('Cell =',num2str(i),' removed with length = ',num2str(li(xvec,i)));
                        xvec = [xvec(1:i-1) xvec(i+1:N+1)];
                        kvec = [kvec(1:i-1) kvec(i+1:N)];
                        a = [a(1:i-1) a(i+1:N)];
                        dead = [dead(1:i-1) dead(i+1:N)];
                        N  = N-1; %Reduce number of cells.

                        if multiple_cell
                        %Adjust differing cell population distribution
                            if i<= border %If cell death occurs to left of border, shift
                                          %border coord left, and reassign populations to
                                          %correct space.
                                border=border-1; 
                                
                            end
                        end
                    end
                    i = i+1; %Step cell-border index.
                end

                %Stochastic process - Birth
                [Lt,L_index] = L_t(xvec,l_birth); %Number of and vector of indices of 
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
                    strcat("Birth at index = ",num2str(i)); 
                    %Split cell and create two copies.
                    midpoint =  (xvec(i) + xvec(i+1))/2;
                    %Introduce first copy at left x-value, and second copy at
                    %midpoint x-value. Reassign x-values.
                    xvec =  [xvec(1:i) midpoint xvec(i+1:N+1)]; %this is now N+2 long
                    kvec = [kvec(1:i) kvec(i) kvec(i+1:N)]; %k_i and a_i are copied.
                    a = [a(1:i) a(i) a(i+1:N)];
                    dead = [dead(1:i) 0 dead(i+1:N)];
                    N = N+1; %Increase number of cells count
                    if multiple_cell %If cell birth occurs to left of border, shift
                                     %border coord right, and reassign populations to
                                       %correct space
                        if i<= border
                            border=border+1;
                        end
                    end
                    %If cell birth occurs to right of border, adjusting number of
                    %cells reassigns populations to correct space.
                end

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
                    a(i) = 0; %Set cell length to zero.
                    kvec(i) = k_death*kvec(i); %Increase spring constant
                    dead(i) = 1;
                    strcat('Cell death occurs for cell=',num2str(i));
                end

                %Account of density after stochastic processes
                if multiple_cell
                    density(1,j+1) = border;
                    density(2,j+1) = N-border;
                else
                    density(j+1) = N;
                end

                %Account of dead cells
                if multiple_cell
                    death_count_1(j) = sum(dead(1:border));
                    death_count_2(j) = sum(dead(border+1:end));
                else
                    death_count_1(j) = sum(dead(:));
                end
            end
        end
        break_val=1;
    end
    index1 = find((density(2,:)==0));
    index2 = find((density(1,:)==0));
    if isempty(index1) == 0
        time_index = index1(1);
        final_time1 = time_index*dt - 100;
    else
        final_time1 = 0;
    end
    if isempty(index2) == 0
        time_index = index2(1);
        final_time2 = time_index*dt - 100;
    else
        final_time2 = 0;
    end
    
end




