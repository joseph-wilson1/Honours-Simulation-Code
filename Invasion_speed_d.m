%% Invasion Speed (as a function of parameters) - d
'Start'
clear all;

% Set up parameters for invasion.
k1 = 20; N = 50; L = 10;
d = 0.00001:0.00001:0.0002;
n = 1000000; nb = 100000;
b = [0.03,0.07]; k2 = [10,15]; dt = 1e-3;
% d = 0.000000005; b = 0.0000000005; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 1;
speeds = zeros(length(b),tests);
speed_ave = zeros(length(b),1);

qtol = 1e-4;
densityfit_speed_pos = 0;
plot_density_snapshot = 0;

fd = figure('Visible','on');
colors = ['r' 'b'; 'g' 'k'];
for k_i = 1:length(k2)
    for b_i = 1:length(b)
        % Run test on invasion speed
        for i = 1:length(d)
            for j = 1:tests
        
                % If gradual death
                if ~id
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function(k1,k2(k_i),L,N,n,nb,b(b_i),d(i),dt,1);
        
                % If instant death   
                else
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function_id(k1,k2(k_i),L,N,n,nb,b(b_i),d(i),dt,1);
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
                    s = strcat("No winner, test ",num2str(j),"d = ",num2str(d(i)));
                    disp(s)
                    speeds(i,j) = 0;
                end

                perc = 100*((i-1)*tests + j)/(length(d)*tests);
                disp(strcat(num2str(round(perc,1)),"% finished"));
            end
        
            % Average data over all tests.
            speed_ave(i) = mean(speeds(i,:));
        end
            
        figure(fd);
        hold on
        plot(d,speed_ave,'+-')
        xlabel("$d$",Interpreter="latex")
        ylabel("$v$",Interpreter="latex")
        ylim([0 inf])
    end
end
legend(strcat("$k_2=$",num2str(k2(1)),", $b=$",num2str(b(1))),strcat("$k_2=$",...
    num2str(k2(1)),", $b=$",num2str(b(2))),strcat("$k_2=$",num2str(k2(2)),...
    ", $b=$",num2str(b(1))),strcat("$k_2=$",num2str(k2(2)),", $b=$",num2str(b(2))),...
    'Interpreter','latex')