function [final_time1,final_time2,SS_density_2,winner,SS_density_1,SS_death,variance_length] = current_invasion_function(k1,k2,L,N,n,nb,b,d,dt,multiple_cell)
%%% k1,k2 are spring constants for pop1 and pop2 resp., L is right spatial boundary, 
%%% N is number of cells, n is number of sims, nb is numbers of sims with hard border,
%%% b is birth prob, d is death prob, dt is time step.

x0 = 0;%N is number of cells.  
border = round(N/2)+1; %Border is an index.
xvec = linspace(x0,L,N+1); %Cells spaced evenly.
kvec = zeros(1,N);
avec = zeros(1,N);

a0 = L/N;
kvec(1:border-1) = k1;
kvec(border:N) = k2;
avec(:) = a0;
dead = zeros(1,N); %Vector to identify wether a cell is dead or not (1,0).

density = zeros(2,n+1);
death_count = zeros(2,n+1);
var_l = zeros(n+1,1);

burn = nb;

winner=0;
final_time1=0;
final_time2=0;

if multiple_cell
    density(1,1) = border-1;
    density(2,1) = N+1-border;
    death_count(1,1) = sum(dead(1:border-1));
    death_count(2,1) = sum(dead(border:end));
    X1 = xvec(1:border);
    X2 = xvec(border:end);
    K1 = kvec(1:border-1);
    K2 = kvec(border:end);
    A1 = avec(1:border-1);
    A2 = avec(border:end);
    D1 = dead(1:border-1);
    D2 = dead(border:end);
else
    nb=0;
    density(1,1) = N;
    death_count(1,1) = sum(dead);
end

plot_on = 1;

eta = 1; 
length_tol = a0/10; % Length at which cells are removed.
lt = 1;
birth_length_tol = lt*a0;

j=0;
t = 0;
k_death = 10; % spring constant for dead cell.

if plot_on
    figure(1)
    clf;
    plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,j,dt)
end

for j=1:nb
    t=t+dt;
    X1_new = hookfun2(X1,K1,A1,dt,eta);
    X1_new = border_test(X1_new);
    [X1_new,K1,A1,D1,~,border] = length_test(X1_new,K1,A1,D1,length_tol,border);
    X1 = X1_new;

    [X1,K1,A1,D1,N1,border] = stochastic_birth(X1,birth_length_tol,b,K1,A1,D1,0,border,t);
    [K1,A1,D1] = stochastic_death(N1,d,K1,A1,D1,k_death);

    X2_new = hookfun2(X2,K2,A2,dt,eta);
    X2_new = border_test(X2_new);
    [X2_new,K2,A2,D2,~,border] = length_test(X2_new,K2,A2,D2,length_tol,border);
    X2 = X2_new;

    [X2,K2,A2,D2,N2,border] = stochastic_birth(X2,birth_length_tol,b,K2,A2,D2,0,border,t);
    [K2,A2,D2] = stochastic_death(N2,d,K2,A2,D2,k_death);

    if multiple_cell
        density(1,j+1) = N1;
        density(2,j+1) = N2;
        death_count(1,j+1) = sum(D1);
        death_count(2,j+1) = sum(D2);
    else
        density(1,j+1) = N;
        death_count(1,j+1) = sum(dead);
        border=N/2;
    end

    if plot_on
        plot_cell([X1(1:end-1) X2],N1+1,k1,k2,x0,L,1,j,dt)
    end
end
if multiple_cell
    xvec = [X1(1:end-1) X2];
    kvec = [K1 K2];
    avec = [A1 A2];
    dead = [D1 D2];
    border = N1+1;
    N = N1+N2;
end

for j=nb:n
    t = t+dt;
    x_old = xvec;
    xvec_new = hookfun2(x_old,kvec,avec,dt,eta);
    xvec_new = border_test(xvec_new);
    [xvec_new,kvec,avec,dead,~,border] = length_test(xvec_new,kvec,avec,dead,length_tol,border);
    xvec=xvec_new;

    if abs(border-1)<1e-4
        final_time1 = (j-nb)*dt;
        final_time2 = 0;
        winner=1;
        break
    elseif abs(border-(N+1))<1e-4
        final_time1 = 0;
        final_time2 = (j-nb)*dt;
        winner=1;
        break
    end

    [xvec,kvec,avec,dead,N,border] = stochastic_birth(xvec,birth_length_tol,b,kvec,avec,dead,multiple_cell,border,t);
    [kvec,avec,dead] = stochastic_death(N,d,kvec,avec,dead,k_death);
    
    if multiple_cell
        density(1,j+1) = border-1;
        density(2,j+1) = N+1-border;
        death_count(1,j+1) = sum(dead(1:border-1));
        death_count(2,j+1) = sum(dead(border:end));
    else
        density(1,j+1) = N;
        death_count(1,j+1) = sum(dead);
        border=N/2;
    end
    
    length_vec = zeros(N,1);
    for i=1:N
        length_vec(i) = li(xvec,i);
    end
    var_l(j+1) = var(length_vec);
    
    if plot_on
        plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,j,dt)
    end
end

%SS_density_1 = density(1,:);
variance_length = mean(var_l(burn:end));
SS_density_1 = mean(density(1,burn:end));
SS_density_2 = mean(density(2,burn:end));
SS_death = mean(death_count(:,burn:end),'all');



















