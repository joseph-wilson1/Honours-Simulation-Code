function [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
    death_count_ss,variance_length,average_density, instant_speed,border_v]...
    = current_invasion_function_id(k1,k2,L,N,n,nb,b,d,dt,multiple_cell)
%%% Inputs:
% k1,k2 are spring constants for pop1 and pop2 resp., L is right spatial boundary, 
% N is number of cells, n is number of sims, nb is numbers of sims with hard border,
% b is birth prob, d is death prob, dt is time step.

%%% Outpus:
% final_time1: Time that population 1 wins, final_time2: time that
% population 2 wins, winner: flag of which population wins, 1 for pop 1, 2
% for pop 2, number_ss_1: average number of cells in population 1 during 
% steady state (after burn), number_ss_2: average number of cells in 
% population 2 during steady state (after burn), death_count_ss: average 
% number of dead cells in simulation during steady state, variance_length:
% average variance in length of cells in simulation, 
% average_density: 2x(n+1) vector, cell i,j contains average density of 
% population i at time-step j. avg density calculated N(t)_i/L_i, 
% instant_speed: 1xn vector, current speed at each time step t,
% calculated as v(t) = (x(border,t)-x(border,t-dt))/dt, 
% border_v: 1x(n-nb+1) vector, border position at each time-step of
% invasion.

% Set up parameters for simulation.
x0 = 0;
border = round(N/2)+1; 
xvec = linspace(x0,L,N+1); 
kvec = zeros(1,N);
avec = zeros(1,N);
a0 = L/N;
kvec(1:border-1) = k1;
kvec(border:N) = k2;
avec(:) = a0; dead = zeros(1,N);

% Set up empty vectors to store information.
number = zeros(2,n+1);
death_count = zeros(2,n+1);
var_l = zeros(n+1,1);
average_density = zeros(2,n+1);
instant_speed = zeros(1,n-nb);
border_v = zeros(1,n-nb+1);

% Density burn number.
burn = nb;

% Set up information variables. 
winner=0;
final_time1=0;
final_time2=0;

% If two populations, set up two initial populations, 
% store initial information.
if multiple_cell
    number(1,1) = border-1;
    number(2,1) = N+1-border;
    average_density(1,1) = number(1,1)/xvec(border);
    average_density(2,1) = number(2,1)/(L-xvec(border));
    death_count(1,1) = sum(dead(1:border-1));
    death_count(2,1) = sum(dead(border:end));
    wall = xvec(border);
    X1 = xvec(1:border);
    X2 = xvec(border:end);
    K1 = kvec(1:border-1);
    K2 = kvec(border:end);
    A1 = avec(1:border-1);
    A2 = avec(border:end);
    D1 = dead(1:border-1);
    D2 = dead(border:end);

% If one population, no wall, store initial information.
else
    nb=0;
    number(1,1) = N;
    death_count(1,1) = sum(dead);
    average_density(1,1) = number(1,1)/L;
end

% Toggle for plotting (try and make this an input).
plot_on = 0;

% Standard simulation parameters (try and make this a txt file).
eta = 1; 
lt = 1;
birth_length_tol = lt*a0;
tol=1e-4;

% Initial time.
t = 0;

% Initial plot.
if plot_on
    figure(1)
    clf;
    plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,0,dt)
end

%%% Order of simulation step (monte carlo step?):
% Time step -> Movement -> Position check -> Stochastic birth 
% -> Stochastic death

%%% If two populations, initial simulation with wall.
for j=1:nb

    % Time-step
    t=t+dt;

    % First population
    X1_new = hookfun2(X1,K1,A1,dt,eta);
    X1 = border_test(X1_new);
    [X1,K1,A1,D1,N1,border] = stochastic_birth(X1,birth_length_tol,b,K1,A1,D1,0,border,t);
    [X1,K1,A1,N1,border] = stochastic_death2(N1,d,X1,K1,A1,border);

    % Second population
    X2_new = hookfun2(X2,K2,A2,dt,eta);
    X2 = border_test(X2_new);
    [X2,K2,A2,D2,N2,border] = stochastic_birth(X2,birth_length_tol,b,K2,A2,D2,0,border,t);
    [X2,K2,A2,N2,border] = stochastic_death2(N2,d,X2,K2,A2,border);

    % Record information at end of time step.
    if multiple_cell
        number(1,j+1) = N1;
        number(2,j+1) = N2;
        death_count(1,j+1) = sum(D1);
        death_count(2,j+1) = sum(D2);
        average_density(1,j+1) = number(1,j+1)/wall;
        average_density(2,j+1) = number(2,j+1)/(L-wall);
    else
        number(1,j+1) = N;
        death_count(1,j+1) = sum(dead);
        border=N/2;
        average_density(1,j+1) = number(1,j+1)/L;
    end
    
    % Plot (if toggle on).
    if plot_on
        plot_cell([X1(1:end-1) X2],N1+1,k1,k2,x0,L,1,j,dt)
    end
end

% If two populations, concatenate into single heterogeneous arrays.
if multiple_cell
    xvec = [X1(1:end-1) X2];
    kvec = [K1 K2];
    avec = [A1 A2];
    dead = [D1 D2];
    border = N1+1;
    N = N1+N2;
    old_border_loc = xvec(border);
    border_v(1,1) = xvec(border);
end

% Simulation without wall.

for j=nb+1:n
    t = t+dt; % Time-step
    xvec_new = hookfun2(xvec,kvec,avec,dt,eta); % Movement
    xvec = border_test(xvec_new); % Position check

    [xvec,kvec,avec,dead,N,border] = stochastic_birth(xvec,birth_length_tol,b,kvec,avec,dead,multiple_cell,border,t); % Stochastic birth
    [xvec,kvec,avec,N,border] = stochastic_death2(N,d,xvec,kvec,avec,border); % Stochastic death

    % Decide if a population has 'won'. Record time of winning invasion.
    if abs(border-1) < tol
        final_time1 = (j-nb)*dt;
        final_time2 = 0;
        winner=1;
        instant_speed(1,j-nb) = (xvec(border)-old_border_loc)/dt;
        border_v(1,(j+1)-nb) = xvec(border);
        break
    elseif abs(border-(N+1)) < tol
        final_time1 = 0;
        final_time2 = (j-nb)*dt;
        winner=2;
        instant_speed(1,j-nb) = (xvec(border)-old_border_loc)/dt;
        border_v(1,(j+1)-nb) = xvec(border);
        break
    end
    
    % Record information at end of time-step.
    if multiple_cell
        number(1,j+1) = border-1;
        number(2,j+1) = N+1-border;
        death_count(1,j+1) = sum(dead(1:border-1));
        death_count(2,j+1) = sum(dead(border:end));
        average_density(1,j+1) = number(1,j+1)/xvec(border);
        average_density(2,j+1) = number(2,j+1)/(L-xvec(border));
        instant_speed(1,j-nb) = (xvec(border)-old_border_loc)/dt;
        border_v(1,(j+1)-nb) = xvec(border);
    else
        number(1,j+1) = N;
        death_count(1,j+1) = sum(dead);
        border=N/2;
        average_density(1,j+1) = number(1,j+1)/L;
    end
    
    % Calculate variation in length.
    length_vec = zeros(N,1);
    for i=1:N
        length_vec(i) = li(xvec,i);
    end
    var_l(j+1) = var(length_vec);

    % Plot (if toggle on).
    if plot_on
        plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,j,dt)
    end

    % Set old border location
    old_border_loc = xvec(border);
end

% End of simulation.

% Calculate outputs.
variance_length = mean(var_l(burn:end));
number_ss_1 = mean(number(1,burn:end));
number_ss_2 = mean(number(2,burn:end));
death_count_ss = mean(death_count(:,burn:end),'all');



















