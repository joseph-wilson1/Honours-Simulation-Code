%% Invasion Speed
'Start'
clear all;

% Set up parameters for invasion.
k1 = 5; N = 50; L = 10;
k2 = 10;
n = 1000000; nb = 100000;
d = 0.0002; b = 0.07; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 1;
speeds = zeros(length(k2),tests);
death_count = zeros(length(k2),tests);
speed_ave = zeros(length(k2),1);
death_ave = zeros(length(k2),1);
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
                instant_speed,border_v,position_array] = ...
                current_invasion_function(k1,k2(i),L,N,n,nb,b,d,dt,1);

        % If instant death   
        else
            [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                death_count_ss,variance_length,average_density,...
                instant_speed,border_v,position_array] = ...
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
    
    if densityfit_speed_pos

        % Linear fit of density vs time
        [~,j] = find(average_density);
        endj = j(end);
        c1 = polyfit((1:endj)*dt,average_density(1,1:endj),1);
        p1 = polyval(c1,(1:endj)*dt);
        c2 = polyfit((1:endj)*dt,average_density(2,1:endj),1);
        p2 = polyval(c2,(1:endj)*dt);
    
        % Plot density fit vs time.
        figure(1)
        hold off
        plot((1:endj)*dt,p1)
        hold on
        plot((1:endj)*dt,p2)
        xlim([dt endj*dt])
        ylim([6 9])
        legend(strcat("Pop. 1, $k_1~ = ~$",num2str(k1)," $~k_2 ~= ~$",num2str(k2(i))),strcat("Pop. 2, $k_1~ = ~$",num2str(k1),", $k_2 ~= ~$",num2str(k2(i))),'Interpreter','latex')

        % Plot instant invasion speed
        figure(2)
        speed_endj = endj-nb;
        plot(((nb+1):endj)*dt,instant_speed(1:speed_endj))

        % Plot border position
        figure(3)
        plot(((nb):endj)*dt,border_v(1:speed_endj+1))
        c3 = polyfit((nb:endj)*dt,border_v(1:speed_endj+1),1);
        p3 = polyval(c3,(nb:endj)*dt);
        hold on
        plot((nb:endj)*dt,p3)
        xlim([nb*dt endj*dt])
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

%% Approx number of cell vs Empirical Population 1
f = find(~cellfun(@isempty,position_array));
f_end = f(end);
position_array = position_array(1:f_end);
clf;
f5=figure('Visible','on');
figure(f5);
N1_em = zeros(1,f_end-nb);
N1_th = zeros(1,f_end-nb);
for l=nb+1:1:f_end
    xvec = position_array{l};
    border_ind = xvec(end);
    N1_em(1,l-nb) = border_ind-1;
    N1_th(1,l-nb) = N_approx(L/2,N/2,xvec(border_ind),k1,b/dt);
end


% c51 = polyfit((nb+1:f_end)*dt,N1_em,1);
% p51 = polyval(c51,(nb+1:f_end)*dt);
% c52 = polyfit((nb+1:f_end)*dt,N1_th,1);
% p52 = polyval(c52,(nb+1:f_end)*dt);
% 
% plot((nb+1:f_end)*dt,p51)
% hold on
% plot((nb+1:f_end)*dt,p52)

plot((nb+1:f_end)*dt,N1_em)
hold on
plot((nb+1:f_end)*dt,N1_th)

%% Approx number of cell vs Empirical Population 2
f = find(~cellfun(@isempty,position_array));
f_end = f(end);
position_array = position_array(1:f_end);
clf;
f5=figure('Visible','on');
figure(f5);
N2_em = zeros(1,f_end-nb);
N2_th = zeros(1,f_end-nb);
for l=nb+1:1:f_end
    xvec = position_array{l};
    border_ind = xvec(end);
    N2_em(1,l-nb) = length(xvec)-border_ind;
    N2_th(1,l-nb) = N_approx(L/2,N/2,L-xvec(border_ind),k2,b/dt);
end


% c61 = polyfit((nb+1:f_end)*dt,N1_em,1);
% p61 = polyval(c61,(nb+1:f_end)*dt);
% c62 = polyfit((nb+1:f_end)*dt,N1_th,1);
% p62 = polyval(c62,(nb+1:f_end)*dt);
% 
% plot((nb+1:f_end)*dt,p61)
% hold on
% plot((nb+1:f_end)*dt,p62)

plot((nb+1:f_end)*dt,N2_em)
hold on
plot((nb+1:f_end)*dt,N2_th)
%% Number of Cell Plot
f = find(~cellfun(@isempty,position_array));
f_end = f(end);
position_array = position_array(1:f_end);
clf;
f5=figure('Visible','on');
figure(f5);
N = zeros(1,f_end-nb);
for l=nb+1:1:f_end
    xvec = position_array{l};
    border_ind = xvec(end);
    xvec = xvec(1:end-1);
    N(1,l-nb) = length(xvec)-1;
end
c5 = polyfit((nb+1:f_end)*dt,N,1);
p5 = polyval(c5,(nb+1:f_end)*dt);
plot((nb+1:f_end)*dt,N)
hold on
plot((nb+1:f_end)*dt,p5)

%% Density profile video
figure(5)
f = find(~cellfun(@isempty,position_array));
f_end = f(end);
position_array = position_array(1:f_end);
for l=nb:1:(nb+100)
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
end

















