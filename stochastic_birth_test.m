'Start'
x0 = 0; L = 10; N = 10; %N is number of cells.  
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

multiple_cell = 1;

b = 0.7;
d = 0.005;

eta = 1; 
n=2; % Number of sims. 
dt = 1e-2; % Time step.
length_tol = a0/10; % Length at which cells are removed.
lt = 1.2;
birth_length_tol = lt*a0;

figure(1)
clf;

k_death = 10; % Factor increase of spring constant for dead cell.

hold on
Ys = zeros(length(xvec(1:border)),1);
plot(xvec(1:border),Ys,'b--o')
title(strcat("Length-dependent birth mechanism - k1 = ",num2str(k1)," k2 = ",num2str(k2)))
xlabel("x")
ylabel("time")
xlim([x0 L])
%Replot second population of cells
Yd = zeros(length(xvec(border:N+1)),1);
plot(xvec(border:N+1),Yd,'r--o')
xlim([x0 L])
drawnow

%%
for j=1:n
    x_old = xvec;
    xvec_new = hookfun2(x_old,kvec,avec,dt,eta);
    xvec_new = border_test(xvec_new);
    [xvec_new,kvec,avec,dead,~,border] = length_test(xvec_new,kvec,avec,dead,length_tol,border);
    xvec=xvec_new;

    [xvec,kvec,avec,dead,N,border] = stochastic_birth(xvec,birth_length_tol,b,kvec,avec,dead,multiple_cell,border);
    %[N,kvec,avec,dead] = stochastic_death(N,d,kvec,avec,dead,k_death);
    
    hold on
    Ys = zeros(length(xvec(1:border)),1)+j*dt;
    plot(xvec(1:border),Ys,'b--o')
    title(strcat("Length-dependent birth mechanism - k1 = ",num2str(1)," k2 = ",num2str(1)))
    xlabel("x")
    ylabel("time")
    xlim([x0 L])
    %Replot second population of cells
    Yd = zeros(length(xvec(border:N+1)),1)+j*dt;
    plot(xvec(border:N+1),Yd,'r--o')
    xlim([x0 L])
    drawnow
end
