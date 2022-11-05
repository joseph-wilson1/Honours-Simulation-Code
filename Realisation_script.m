%% Script to produce figures.
'Start'
clear all;

% Set up parameters for simulation.
k1 = 10; N = 50; L = 10; k2 = 10;
n = 10000; nb = 1;
d = 0.00005; b = 0.07; dt = 1e-3;
id = 0;

number_ave1 = zeros(1,n+1); % Number of cells vs time (1)

mc=0; % Toggle for single vs multi-population.

% Toggle on plotting flag in invasion_functions.

% If gradual death
if ~id
    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
        death_count_ss,variance_length,average_density,...
        instant_speed,border_v,position_array,number_cells] = ...
        current_invasion_function(k1,k2,L,N,n,nb,b,d,dt,mc);

% If instant death   
else
    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
        death_count_ss,variance_length,average_density,...
        instant_speed,border_v,position_array,number_cells] = ...
        current_invasion_function_id(k1,k2,L,N,n,nb,b,d,dt,mc);
end
number1 = number_cells(1,:);
mean_number = mean(number1(1,round(n/2):n));
%%
figure(3)
hold off
plot((0:n)*dt,number1)
hold on
yline(mean_number,'r','LineWidth',2)
xlim([0 n*dt])
xlabel("$t$",Interpreter="latex")
ylabel("$N(t)$",Interpreter="latex")

