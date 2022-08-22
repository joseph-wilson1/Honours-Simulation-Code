%% Invasion Speed
'Start'
clear all;

% Set up parameters for invasion.
k1 = 5; N = 50; L = 10;
k2 = [5 7 10];
n = 100000; nb = 10000;
d = 0.002; b = 0.7; dt = 1e-2;
id = 0;

% Set up parameters for testing.

tests = 1;
speeds = zeros(length(k2),tests);
death_count = zeros(length(k2),tests);
speed_ave = zeros(length(k2),1);
death_ave = zeros(length(k2),1);
qtol = 1e-4;
densityfit = 1;

% Run test on invasion speed
for i = 1:length(k2)
    for j = 1:tests

        % If gradual death
        if ~id
            [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                death_count_ss,variance_length,average_density] = ...
                current_invasion_function(k1,k2(i),L,N,n,nb,b,d,dt,1);

        % If instant death   
        else
            [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                death_count_ss,variance_length,average_density] = ...
                current_invasion_function_id(k1,k2(i),L,N,n,nb,b,d,dt,1);
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
            s = strcat("No winner, test ",num2str(j),"k_2 = ",num2str(k2(i)));
            disp(s)
            speeds(i,j) = 0;
        end

    end

    % Average speed over all tests.
    speed_ave(i) = mean(speeds(i,:));
    
    if densityfit

        % Linear fit of density vs time
        [~,j] = find(average_density);
        endj = j(end);
        c1 = polyfit(1:endj,average_density(1,1:endj),1);
        p1 = polyval(c1,1:endj);
        c2 = polyfit(1:endj,average_density(2,1:endj),1);
        p2 = polyval(c2,1:endj);
    
        % Plot density fit vs time.
        figure(i)
        hold off
        plot((1:endj)*dt,p1)
        hold on
        plot((1:endj)*dt,p2)
        xlim([dt endj*dt])
        ylim([6 8])
        legend(strcat("Pop. 1, $k_1~ = ~$",num2str(k1)," $~k_2 ~= ~$",num2str(k2(i))),strcat("Pop. 2, $k_1~ = ~$",num2str(k1),", $k_2 ~= ~$",num2str(k2(i))),'Interpreter','latex')
    end
end