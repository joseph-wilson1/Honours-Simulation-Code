%%% Trying out a border system

'Start'
x0 = 0; L = 10; N = 50; %N is number of cells.  
border = round(N/2)+1; %Border is an index.
xvec = linspace(x0,L,N+1); %Cells spaced evenly.
kvec = zeros(1,N);
avec = zeros(1,N);

k1 = 5; k2 = 5;

a0 = L/N;
kvec(1:border-1) = k1;
kvec(border:N) = k2;
avec(:) = a0;
dead = zeros(1,N); %Vector to identify wether a cell is dead or not (1,0).

n=100000; % Number of sims. 
nb = 100; % Number of sims with hard border.
dt = 1e-2; % Time step.

density = zeros(2,n+1);
death_count = zeros(2,n+1);

multiple_cell = 0;

plot_on = 0;

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

b = 0.7;
d = 0.001;

eta = 1; 
length_tol = a0/10; % Length at which cells are removed.
lt = 1.2;
birth_length_tol = lt*a0;

figure(1)
clf;
t = 0;
j=0;
k_death = 2; % Factor increase of spring constant for dead cell.

if plot_on
    figure(1)
    clf;
    plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,j,dt)
end
%%
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
        plot_cell([X1(1:end-1) X2],N1+1,k1,k2,x0,L,multiple_cell,j,dt)
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

    [xvec,kvec,avec,dead,N,border] = stochastic_birth(xvec,birth_length_tol,b,kvec,avec,dead,multiple_cell,border,t);
    [kvec,avec,dead] = stochastic_death(N,d,kvec,avec,dead,k_death);
    
    if border==1
        break
    elseif border==N+1
        break
    end

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
    
    if plot_on
        plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,j,dt)
    end
end

%%
plot_density(density,n,dt,multiple_cell,k1,k2)
burn = n/2;
SS_density_1 = mean(density(1,burn:end));
SS_density_2 = mean(density(2,burn:end));
plot_dead(death_count,n,dt,k1,k2,multiple_cell)
%final
figure(1)
plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,j,dt)





















