%%% 02/05/22 - Current issue is bias towards first population. I have no
%%% idea why.

'Start'
x0 = 0; L = 10; N = 50; %N is number of cells.  
border = round(N/2)+1; %Border is an index.
xvec = linspace(x0,L,N+1); %Cells spaced evenly.
kvec = zeros(1,N);
avec = zeros(1,N);

k1 = 5; k2 = .5;

a0 = L/N;
kvec(1:border-1) = k1;
kvec(border:N) = k2;
avec(:) = a0;
dead = zeros(1,N); %Vector to identify wether a cell is dead or not (1,0).

n=100000; % Number of sims. 
dt = 1e-2; % Time step.

density = zeros(2,n+1);
death_count = zeros(2,n+1);

multiple_cell = 1;

if multiple_cell
    density(1,1) = border-1;
    density(2,1) = N+1-border;
    death_count(1,1) = sum(dead(1:border-1));
    death_count(2,1) = sum(dead(border:end));
else
    density(1,1) = N;
    death_count(1,1) = sum(dead);
end

b = 0.7;
d = 0.001;

plot_on = 0;

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
for j=1:n
    t = t+dt;
    x_old = xvec;
    xvec_new = hookfun2(x_old,kvec,avec,dt,eta);
    xvec_new = border_test(xvec_new);
    [xvec_new,kvec,avec,dead,~,border] = length_test(xvec_new,kvec,avec,dead,length_tol,border);
    xvec=xvec_new;

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
    end
    
    if plot_on
        plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,j,dt)
    end
end

%%
plot_density(density,n,dt,multiple_cell,k1,k2)
burn = 100;
SS_density_1 = mean(density(1,burn:end));
SS_density_2 = mean(density(2,burn:end));
plot_dead(death_count,n,dt,k1,k2,multiple_cell)
%final
plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,j,dt)






















