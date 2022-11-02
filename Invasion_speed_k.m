%% Invasion Speed (as a function of parameters) - k1
'Start'
clear all;

% Set up parameters for invasion.
k1 = 20; N = 50; L = 10;
k2 = 5:0.5:20;
n = 1000000; nb = 100000;
d = [0.00005,0.0002]; b = [0.03,0.07]; dt = 1e-3;
% d = 0.000000005; b = 0.0000000005; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 1;
speeds = zeros(length(k2),tests);
death_count = zeros(length(k2),tests);
speed_ave = zeros(length(k2),1);
death_ave = zeros(length(k2),1);
border_ave = zeros(tests,n-nb+1);
density_ave1 = zeros(tests,n+1);
density_ave2 = zeros(tests,n+1);
number_ave1 = zeros(tests,n+1);
number_ave2 = zeros(tests,n+1);
end_time = zeros(tests,1);

qtol = 1e-4;
densityfit_speed_pos = 0;
plot_density_snapshot = 0;

f9 = figure('Visible','on');
colors = ['r' 'b'; 'g' 'k'];
for b_i = 1:length(b)
    for d_i = 1:length(d)
        % Run test on invasion speed
        for i = 1:length(k2)
            for j = 1:tests
        
                % If gradual death
                if ~id
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function(k1,k2(i),L,N,n,nb,b(b_i),d(d_i),dt,1);
        
                % If instant death   
                else
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function_id(k1,k2(i),L,N,n,nb,b(b_i),d(d_i),dt,1);
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
                
                % Record data from test j
                border_ave(j,:) = border_v;
                density_ave1(j,:) = average_density(1,:);
                density_ave2(j,:) = average_density(2,:);
                number_ave1(j,:) = number_cells(1,:);
                number_ave2(j,:) = number_cells(2,:);
                et = find(average_density(1,:));
                end_time(j,1) = et(end);
                
                perc = 100*((i-1)*tests + j)/(length(k2)*tests);
                disp(strcat(num2str(round(perc,1)),"% finished"));
            end
        
            % Average data over all tests.
            speed_ave(i) = mean(speeds(i,:));
            border_ave = mean(border_ave,1);
            density_ave1 = mean(density_ave1,1);
            density_ave2 = mean(density_ave2,1);
            number_ave1 = mean(number_ave1,1);
            number_ave2 = mean(number_ave2,1);
            
            if densityfit_speed_pos
        %%
                % Plot average density vs time.
                [~,j1] = find(density_ave1);
                endj1 = j1(end);
                [~,j2] = find(density_ave2);
                endj2 = j2(end);
                figure(1)
                hold off
                plot((1:endj1)*dt,density_ave1(1,1:endj1))
                hold on
                plot((1:endj2)*dt,density_ave2(1,1:endj2))
                xlim([dt max([endj1,endj2])*dt])
                ylim([0 20])
                legend(strcat("Pop. 1, $k_1~ = ~$",num2str(k1)," $~k_2 ~= ~$",num2str(k2(i))),strcat("Pop. 2, $k_1~ = ~$",num2str(k1),", $k_2 ~= ~$",num2str(k2(i))),'Interpreter','latex')
        %%
                % Plot border position
                figure(2)
                speed_endj = endj1-nb;
                hold off
                plot(((nb):endj1-1)*dt,border_ave(1:speed_endj),'LineWidth',3)
                c3 = polyfit((nb:round(speed_endj/2,0)+nb)*dt,border_ave(1:round(speed_endj/2,0)+1),1);
                p3 = polyval(c3,(nb:endj1)*dt);
                hold on
                plot((nb:endj1)*dt,p3)
                xlabel("Time $t$",Interpreter="latex")
                ylabel("Border Position $x_i(t)$",Interpreter="latex")
                xlim([nb*dt endj1*dt])
                ylim([0 5])
                title(strcat("Gradient of Linear Curve = ",num2str(round(c3(1),3))))
            end
        
        end
        figure(f9);
        hold on
        plot(k2,speed_ave,'o-')
        xlabel("$k_2$",Interpreter="latex")
        ylabel("$v$",Interpreter="latex")
        ylim([0 inf])
    end
end
legend(strcat("$b=$",num2str(b(1)),", $d=$",num2str(d(1))),strcat("$b=$",...
    num2str(b(1)),", $d=$",num2str(d(2))),strcat("$b=$",num2str(b(2)),...
    ", $d=$",num2str(d(1))),strcat("$b=$",num2str(b(2)),", $d=$",num2str(d(2))),...
    'Interpreter','latex')


%% invasion approximation #1
a0=0.2;
v_approx = zeros(1,length(k2));
for i=1:length(k2)
    Ni1 = N_approx(L/2,N/2,L/2,k1,b(1)/dt);
    Ni2 = N_approx(L/2,N/2,L/2,k2(i),b(1)/dt);
    li1 = L/2/Ni1; li2 = L/2/Ni2;
    v_approx(1,i) = -k1*(li1-a0) + k2(i)*(li2-a0);
end
hold on
plot(k2,v_approx)
legend("Numerical","Approximation")















