%% Script to produce movie of invasion process.
'Start'
clear all;

% Set up parameters for invasion.
k1 = 20; N = 50; L = 10;
k2 = 10;
n = 2000000; nb = 200000;
d = 0.00005; b = 0.03; dt = 1e-3;
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

%% Cell position plot
f = find(~cellfun(@isempty,position_array));
f_end = f(end);
position_array = position_array(1:f_end);
clf;
f4=figure('Visible','off');
f4.Position = [210 623 1500 100];
figure(f4);
movieSize = 4000;
m_i = 1;
for l=nb+1:1:nb+movieSize+1
    xvec = position_array{l};
    border_ind = xvec(end);
    xvec = xvec(1:end-1);
    y = zeros(1,length(xvec));
    hold off
    plot(xvec(1:border_ind),y(1:border_ind),'b:|','LineWidth',1)
    hold on
    plot(xvec(border_ind:end),y(border_ind:end),'r:|','LineWidth',1)
    ylim([-0.1 0.1])
    hold on
    plot(xvec(border_ind),y,'k|')
    title(strcat("$t~=~$",num2str(round((l-nb)*dt,2)),", $v_t~=~$",num2str(round(instant_speed((l-nb)),3)),", $N_1~=~$",num2str(border_ind-1),", $N_2~=~$",num2str(length(xvec)-border_ind)),strcat("Border Position $~=~$",num2str(xvec(border_ind))),'Interpreter','latex')
    
    movieVector(m_i) = getframe(f4);
    m_i=m_i+1;
end

myWriter = VideoWriter('cells');
myWriter.FrameRate = 60;

open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);










