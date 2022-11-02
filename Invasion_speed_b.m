%% Invasion Speed (as a function of parameters) - b 
'Start'
clear all;

% Set up parameters for invasion.
k1 = 20; N = 50; L = 10;
b = 0.005:0.005:0.1;
n = 1000000; nb = 100000;
d = [0.00005, 0.0002]; k2 = [10,15]; dt = 1e-3;
% d = 0.000000005; b = 0.0000000005; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 20;
speeds = zeros(length(b),tests);
speed_ave = zeros(length(b),1);

qtol = 1e-4;
densityfit_speed_pos = 0;
plot_density_snapshot = 0;

f9 = figure('Visible','on');
colors = ['r' 'b'; 'g' 'k'];
for k_i = 1:length(k2)
    for d_i = 1:length(d)
        % Run test on invasion speed
        for i = 1:length(b)
            for j = 1:tests
        
                % If gradual death
                if ~id
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function(k1,k2(k_i),L,N,n,nb,b(i),d(d_i),dt,1);
        
                % If instant death   
                else
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function_id(k1,k2(k_i),L,N,n,nb,b(i),d(d_i),dt,1);
                end
        
                % Record invasion speed of winner -> + for pop1, - for pop2.
                if winner
                    if abs(winner-1) < qtol
                        speeds(i,j) = (L/2) / final_time1;
                    elseif abs(winner-2) < qtol
                        speeds(i,j) = -(L/2) / final_time2;
                    end
        
                % If no winner, speed = 0.
                else
                    s = strcat("No winner, test ",num2str(j),"k_2 = ",num2str(k2(k_i)));
                    disp(s)
                    speeds(i,j) = 0;
                end
                
                perc = 100*((i-1)*tests + j)/(length(b)*tests);
                disp(strcat(num2str(round(perc,1)),"% finished"));
            end
        
            % Average data over all tests.
            speed_ave(i) = mean(speeds(i,:));   
        end
        figure(f9);
        hold on
        plot(b,speed_ave,'*-')
        xlabel("$b$",Interpreter="latex")
        ylabel("$v$",Interpreter="latex")
        ylim([0 inf])
    end
end
legend(strcat("$k_2=$",num2str(k2(1)),", $d=$",num2str(d(1))),strcat("$k_2=$",...
    num2str(k2(1)),", $d=$",num2str(d(2))),strcat("$k_2=$",num2str(k2(2)),...
    ", $d=$",num2str(d(1))),strcat("$k_2=$",num2str(k2(2)),", $d=$",num2str(d(2))),...
    'Interpreter','latex')