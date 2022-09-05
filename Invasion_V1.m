%% Invasion Speed
'Start'
clear all;

% Set up parameters for invasion.
k1 = 5; N = 50; L = 10;
k2 = 7;
n = 1000000; nb = 100000;
d = 0.00005; b = 0.07; dt = 1e-3;
% d = 0.000000005; b = 0.0000000005; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 20;
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
densityfit_speed_pos = 1;
plot_density_snapshot = 0;

% Run test on invasion speed
for i = 1:length(k2)
    for j = 1:tests

        % If gradual death
        if ~id
            [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                death_count_ss,variance_length,average_density,...
                instant_speed,border_v,position_array,number_cells] = ...
                current_invasion_function(k1,k2(i),L,N,n,nb,b,d,dt,1);

        % If instant death   
        else
            [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                death_count_ss,variance_length,average_density,...
                instant_speed,border_v,position_array,number_cells] = ...
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
        
        % Record data from test j
        border_ave(j,:) = border_v;
        density_ave1(j,:) = average_density(1,:);
        density_ave2(j,:) = average_density(2,:);
        number_ave1(j,:) = number_cells(1,:);
        number_ave2(j,:) = number_cells(2,:);
        et = find(average_density(1,:));
        end_time(j,1) = et(end);

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

%% Cell position plot
f = find(~cellfun(@isempty,position_array));
f_end = f(end);
position_array = position_array(1:f_end);
clf;
f4=figure('Visible','on');
f4.Position = [210 623 1500 100];
figure(f4);
for l=nb+1:10:f_end
    xvec = position_array{l};
    border_ind = xvec(end);
    xvec = xvec(1:end-1);
    y = zeros(1,length(xvec));
    hold off
    plot(xvec,y,':|')
    ylim([-0.1 0.1])
    hold on
    plot(xvec(border_ind),y,'r|')
%     title(num2str(instant_speed(l-nb)))
    title(strcat("$t~=~$",num2str(l*dt),", $v_t~=~$",num2str(instant_speed((l-nb))),", $N_1~=~$",num2str(border_ind-1),", $N_2~=~$",num2str(length(xvec)-border_ind)," $N_{1,approx} = $",num2str(N_approx(L/2,N/2,xvec(border_ind),k1,b/dt)))...
        ,strcat("Border Position $~=~$",num2str(xvec(border_ind)),", $D_1~=~$",num2str(average_density(1,l)),", $D_2~=~$",num2str(average_density(2,l)),", $D_{r,t}=\frac{D_1}{D_2}~=~$",...
        num2str(average_density(1,l)/average_density(2,l)),", $D_{r,nb} = $",num2str(average_density(1,nb)/average_density(2,nb)),", Number of cells $~=~$", num2str(length(xvec)-1)),'Interpreter','latex')
    drawnow()
end

%% Number of Cells Empirical vs Theoretical
[~,f] = find(number_ave1);
f_end = f(end);
N1_th = zeros(1,f_end-nb);
N2_th = zeros(1,f_end-nb);
for l=nb+1:1:f_end
    N1_th(1,l-nb) = N_approx(L/2,N/2,border_ave(l-nb),k1,b/dt);
    N2_th(1,l-nb) = N_approx(L/2,N/2,L-border_ave(l-nb),k2,b/dt);
end
figure(5)
hold off
plot((nb+1:f_end)*dt,number_ave1(nb+1:f_end),'r')
hold on
plot((nb+1:f_end)*dt,N1_th,'k--')
plot((nb+1:f_end)*dt,number_ave2(nb+1:f_end),'b')
hold on
plot((nb+1:f_end)*dt,N2_th,'k--')
xlabel("Time $t$",Interpreter="latex")
ylabel("$N(t)$",Interpreter="latex")
xlim([(nb)*dt f_end*dt])
legend(strcat("Pop. 1, $k_1~ = ~$",num2str(k1)),'',strcat("Pop. 2, $k_2 ~= ~$",num2str(k2(i))),'','Interpreter','latex')

%% Density profile video

f = find(~cellfun(@isempty,position_array));
f_end = f(end);
position_array = position_array(1:f_end);
clf;
f5=figure('Visible','on');
figure(f5);
for l=nb:10:f_end
    xvec = position_array{l};
    border_ind = xvec(end);
    xvec = xvec(1:end-1);
    ave_dens = rolling_average_density(xvec);
    hold off
    plot(xvec,[ave_dens ave_dens(end)])
    ylim([0 20])
    hold on
    plot(xvec(border_ind),ave_dens(border_ind),'*r')
    legend(strcat("$t~=~$",num2str(l*dt)),strcat('$x_{border}~=~$',num2str(xvec(border_ind))),'Interpreter','latex')
    drawnow()
end

















