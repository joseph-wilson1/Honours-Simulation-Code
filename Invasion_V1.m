%% Invasion Speed
'Start'
clear all;

% Set up parameters for invasion.
k1 = 10; N = 100; L = 20;
k2 = [20];
n = 1000000; nb = 100000;
d = [0.00005]; b = [0.07]; dt = 1e-3;
% d = 0.000000005; b = 0.0000000005; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 1; % Number of tests
speeds = zeros(length(k2),tests); % Invasion speed for k2(i) in test j
death_count = zeros(length(k2),tests); % Not in use.
speed_ave = zeros(length(k2),1); % Average invasion speed for k2(i)
death_ave = zeros(length(k2),1); % Not in use.
border_ave = zeros(tests,n-nb+1); % Border positions vs time in each test.
density_ave1 = zeros(tests,n+1); % Average density vs time in each test (1)
density_ave2 = zeros(tests,n+1); % Average density vs time in each test (2)
number_ave1 = zeros(tests,n+1); % Number of cells vs time in each test (1)
number_ave2 = zeros(tests,n+1); % Number of cells vs time in each test (2)
end_time = zeros(tests,1); % End iteration of invasion in each test.

qtol = 1e-4;
averagedensity_plot = 1;
border_plot = 1;
plot_density_snapshot = 0;
mc=1;

% Run test on invasion speed
for m=1:length(d)
    for p=1:length(b)
        for i = 1:length(k2)
            for j = 1:tests
        
                % If gradual death
                if ~id
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function(k1,k2(i),L,N,n,nb,b(p),d(m),dt,mc);
        
                % If instant death   
                else
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function_id(k1,k2(i),L,N,n,nb,b(p),d(m),dt,mc);
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
            
            if averagedensity_plot
        
                % Plot average density vs time.
                [~,j1] = find(density_ave1);
                endj1 = j1(end);
                [~,j2] = find(density_ave2);
                endj2 = j2(end);
                figure(1)
                hold off
                plot((1:endj1)*dt,density_ave1(1,1:endj1),'LineWidth',2)
                hold on
                plot((1:endj2)*dt,density_ave2(1,1:endj2),'LineWidth',2)
                xlim([dt max([endj1,endj2])*dt])
                ylim([0 20])
                xlabel("Time $t$",Interpreter="latex")
                ylabel("$\rho_{avg,i}(t)$",Interpreter="latex")
                legend(strcat("Pop. $i=1$, $k_1~ = ~$",num2str(k1)),strcat("Pop. $i=2$, $k_2 ~= ~$",num2str(k2(i))),'Interpreter','latex')
            end
%%
            if border_plot
        
                % Plot border position
                figure(i+1)
                speed_endj = endj1-nb+1;
                hold on
                plot(((nb):endj1)*dt,border_ave(1:speed_endj),'LineWidth',3)
                c3 = polyfit((nb:round(speed_endj/3,0)+nb)*dt,border_ave(1:round(speed_endj/3,0)+1),1);
                p3 = polyval(c3,(nb:endj1)*dt);
                hold on
                plot((nb:endj1)*dt,p3)
                xlabel("Time $t$",Interpreter="latex")
                ylabel("Border Position $S(t)$",Interpreter="latex")
%                 xlim([nb*dt endj1*dt])
                ylim([0 L/2])
%                 title(strcat("Gradient of Linear Curve = ",num2str(round(c3(1),3))))
            end
        end
    end
end
%Add legend to figure(2)/figure(3)
figure(2)
legend(strcat("$d~=~$",num2str(d(1)),", $b~=~$",num2str(b(1))),'',...
    strcat("$d~=~$",num2str(d(1)),"$, b~=~$",num2str(b(2))),'', ...
    strcat("$d~=~$",num2str(d(2)),"$, b~=~$",num2str(b(1))),'',...
    strcat("$d~=~$",num2str(d(2)),"$, b~=~$",num2str(b(2))),'',...
    "interpreter","latex");
figure(3)
legend(strcat("$d~=~$",num2str(d(1)),"$, b~=~$",num2str(b(1))),'',...
    strcat("$d~=~$",num2str(d(1)),"$, b~=~$",num2str(b(2))),'', ...
    strcat("$d~=~$",num2str(d(2)),"$, b~=~$",num2str(b(1))),'',...
    strcat("$d~=~$",num2str(d(2)),"$, b~=~$",num2str(b(2))),'',...
    "interpreter","latex");

%% Add vline to figure(2)
figure(2)
hold on
xline(300,'LineWidth',1)
xlim([100 650])
%% Modify density plot
figure(1)
hold on
ylim([0 12])
xline(100,'LineWidth',1)
xline(300,'LineWidth',1)
legend(strcat("Pop. 1, $k_1~ = ~$",num2str(k1)),strcat("Pop. 2, $k_2 ~= ~$",num2str(k2(i))),'','Interpreter','latex')

%% Invasion Speed Plot
hold on
plot(k2,speed_ave,'ro')
xlabel("$k_2$",Interpreter="latex")
ylabel("$v$",Interpreter="latex")
ylim([0 inf])



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
    N1_th(1,l-nb) = N_approx(L/2,N/2,border_ave(l-nb),k1,b(2)/dt);
    N2_th(1,l-nb) = N_approx(L/2,N/2,L-border_ave(l-nb),k2(2),b(2)/dt);
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
density_plotting = 1;
force_plotting = ~density_plotting;
m = 6;
t=tiledlayout(m,1);
% for l=[nb-1000,round(linspace(nb,f_end/2,m-1),0)]
for l=[nb-1000,nb+(1:m-1)*18000]
% for l=nb:10:f_end
    nexttile
    xvec = position_array{l};
    border_ind = xvec(end);
    border_pos = xvec(border_ind);
    xvec = xvec(1:end-1);
    mvec = midp_cell(xvec); % midpoint of cell

    if density_plotting
    %     dvec = rolling_average_density(xvec);
        dvec = density_vector(xvec); % density of cell
        border_dens = (dvec(border_ind-1)+dvec(border_ind))/2;
        hold off
        plot(mvec,dvec,'.')
        hold on
        plot(border_pos,border_dens,'*r')
        ylim([5 10])
        text(17,8.5,strcat("$t~=~$",num2str(round(l*dt,2))),'Interpreter','latex')
    end
    
    if force_plotting
        f = cell_force(xvec,border_ind,k1,k2(end),L/N);
        border_force = (f(border_ind-1)+f(border_ind))/2;
        hold off
        plot(mvec,f,'.r')
        hold on
        plot(border_pos,border_force,'*k')
        ylim([0 2.5])
        text(17,1.8,strcat("$t~=~$",num2str(round(l*dt,2))),'Interpreter','latex')
    end



    set(gca, 'XTick', []);
end
t.TileSpacing = 'none';
set(gca, 'XTick', 0:20)
if density_plotting
    ylabel(t,"$\rho(x,t)$",Interpreter="latex")
elseif force_plotting
    ylabel(t,"$f(x,t)$",Interpreter="latex")
end
xlabel(t,"$x$",Interpreter="latex")
%% Space-time plot (single-cell)
f = find(~cellfun(@isempty,position_array));
f_end = f(end);
position_array = position_array(1:f_end);
clf;
f6=figure('Visible','on');
figure(f6);
t=0;
for l=1:f_end
    xvec = position_array{l};
    border_ind = xvec(end);
    xvec = xvec(1:end-1);
    y = zeros(1,length(xvec))+(l-1)*dt;
    hold on
    plot(xvec,y,'b.')
    drawnow()
end
xlabel('x')
ylabel('t')