%%% Script to produce number of cells at steady state plots.
plot_approx = 1; % Toggle to plot approx.
run_test = 0;
%%
% Number of Cells at Steady State (Single-Cell) - function of k
'Start'

% Set up parameters for simulation.
N = 50; L = 10;
k2 = 0.5:0.5:20;
n = 100000; nb = 10000;
d = [0.00005, 0.0002]; b = [0.03, 0.07]; dt = 1e-3;
id = 0;

% Set up parameters for testing.

tests = 20;
ss_number = zeros(length(k2),1);
ss_number_test = zeros(tests,1);

mc=0;
if run_test
    f1 = figure('Visible','on');
    for b_i = 1:length(b)
        for d_i = 1:length(d)
            % Run test on number
            for i = 1:length(k2)
                for j = 1:tests
            
                    % If gradual death
                    if ~id
                        [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                            death_count_ss,variance_length,average_density,...
                            instant_speed,border_v,position_array,number_cells] = ...
                            current_invasion_function(k2(i),k2(i),L,N,n,nb,b(b_i),d(d_i),dt,mc);
            
                    % If instant death   
                    else
                        [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                            death_count_ss,variance_length,average_density,...
                            instant_speed,border_v,position_array,number_cells] = ...
                            current_invasion_function_id(k2(i),k2(i),L,N,n,nb,b(b_i),d(d_i),dt,mc);
                    end
    
                    % Record data from test j
                    ss_number_test(j,1) = number_ss_1;
                    
                    perc = 100*((i-1)*tests + j)/(length(k2)*tests);
                    disp(strcat(num2str(round(perc,1)),"% finished"));
                end
            
                % Average data over all tests.
                ss_number(i) = mean(ss_number_test,'all');
            
            end
            figure(f1);
            hold on
            plot(k2,ss_number,'o-')
            xlabel("$k$",Interpreter="latex")
            ylabel("Number of Cells",Interpreter="latex")
        end
    end
end
legend(strcat("$b=$",num2str(b(1)),", $d=$",num2str(d(1))),strcat("$b=$",...
    num2str(b(1)),", $d=$",num2str(d(2))),strcat("$b=$",num2str(b(2)),...
    ", $d=$",num2str(d(1))),strcat("$b=$",num2str(b(2)),", $d=$",num2str(d(2))),...
    'Interpreter','latex')


if plot_approx
    % Plot approximation.
    ylim([60 75])
    approx_number = N_approx(L,N,L,k2,b(2)/dt);
    plot(k2,approx_number,'LineWidth',2)
    approx_number = N_approx(L,N,L,k2,b(1)/dt);
    plot(k2,approx_number,'LineWidth',2)
    legend(strcat("$b=$",num2str(b(1)),", $d=$",num2str(d(1))), ...
        strcat("$b=$",num2str(b(1)),", $d=$",num2str(d(2))),...
        strcat("$b=$",num2str(b(2)),", $d=$",num2str(d(1))),...
        strcat("$b=$",num2str(b(2)),", $d=$",num2str(d(2))),...
        "Theoretical Fit b $= 0.07$",...
        "Theoretical Fit b $= 0.03$",...
        'Interpreter','latex')
end

%%
% Number of Cells at Steady State (Single-Cell) - function of d.
'Start'

% Set up parameters for simulation.
N = 50; L = 10;
k2 = [4 10];
n = 100000; nb = 10000;
d = 0.00001:0.00001:0.0002; 
b = [0.03, 0.07]; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 20;
ss_number = zeros(length(d),1);
ss_number_test = zeros(tests,1);

mc=0;

% f2 = figure('Visible','on');
colors = ['r','b','g','k'];
c=0;
for b_i = 1:length(b)
    % Run test on invasion speed
    for k_i = 1:length(k2)        
        c=c+1;
        for i = 1:length(d)
            if run_test
                for j = 1:tests
            
                    % If gradual death
                    if ~id
                        [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                            death_count_ss,variance_length,average_density,...
                            instant_speed,border_v,position_array,number_cells] = ...
                            current_invasion_function(k2(k_i),k2(k_i),L,N,n,nb,b(b_i),d(i),dt,mc);
            
                    % If instant death   
                    else
                        [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                            death_count_ss,variance_length,average_density,...
                            instant_speed,border_v,position_array,number_cells] = ...
                            current_invasion_function_id(k2(k_i),k2(k_i),L,N,n,nb,b(b_i),d(i),dt,mc);
                    end
                    
                    % Record data from test j
                    ss_number_test(j,1) = number_ss_1;
                    
                    perc = 100*((i-1)*tests + j)/(length(d)*tests);
                    disp(strcat(num2str(round(perc,1)),"% finished"));
                end
        
                % Average data over all tests.
                ss_number(i) = mean(ss_number_test,'all');
                
            end
        end

%         figure(f2);
        hold on
%         plot(d,ss_number,'+-','Color',colors(c))
        if plot_approx
            d_N_approx = zeros(1,length(d)) + N_approx(L,N,L,k2(k_i),b(b_i)/dt);
            hold on
            plot(d,d_N_approx,'Color',colors(c))
        end
        xlabel("$d$",Interpreter="latex")
        ylabel("Number of Cells",Interpreter="latex")
    end
end
%%
xlim([d(1) d(end)])
ylim([60 75])
xlabel("$d$",Interpreter="latex")
ylabel("Number of Cells",Interpreter="latex")
legend(strcat("$b=$",num2str(b(1)),", $k_2=$",num2str(k2(1))), ...
    strcat("$b=$",num2str(b(1)),", $k_2=$",num2str(k2(2))), ...
    strcat("$b=$",num2str(b(2)),", $k_2=$",num2str(k2(1))), ...
    strcat("$b=$",num2str(b(2)),", $k_2=$",num2str(k2(2))),'', ...
    '','','',...
    strcat("$b=$",num2str(b(2)),", $k_2=$",num2str(k2(2))),'',...
    'Interpreter','latex')
%%
% Number of Cells at Steady State (Single-Cell) - function of b.
'Start'

% Set up parameters for simulation.
N = 50; L = 10;
k2 = [4 10];
n = 100000; nb = 10000;
d = [0.00005, 0.0002]; 
b = 0.01:0.01:0.1; dt = 1e-3;
id = 0;

% Set up parameters for testing.

tests = 20;
ss_number = zeros(length(b),1);
ss_number_test = zeros(tests,1);

mc=0;
if run_test
    f3 = figure('Visible','on');
    for d_i = 1:length(d)
        % Run test on invasion speed
        for k_i = 1:length(k2)        
            for i = 1:length(b)
                for j = 1:tests
            
                    % If gradual death
                    if ~id
                        [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                            death_count_ss,variance_length,average_density,...
                            instant_speed,border_v,position_array,number_cells] = ...
                            current_invasion_function(k2(k_i),k2(k_i),L,N,n,nb,b(i),d(d_i),dt,mc);
            
                    % If instant death   
                    else
                        [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                            death_count_ss,variance_length,average_density,...
                            instant_speed,border_v,position_array,number_cells] = ...
                            current_invasion_function_id(k2(k_i),k2(k_i),L,N,n,nb,b(i),d(d_i),dt,mc);
                    end
                    
                    % Record data from test j
                    ss_number_test(j,1) = number_ss_1;
                    
                    perc = 100*((i-1)*tests + j)/(length(b)*tests);
                    disp(strcat(num2str(round(perc,1)),"% finished"));
                end
        
                % Average data over all tests.
                ss_number(i) = mean(ss_number_test,'all');
                
            end
    
            figure(f3);
            hold on
            plot(b,ss_number,'*-')
            xlabel("$b$",Interpreter="latex")
            ylabel("Number of Cells",Interpreter="latex")
        end
    end
end
%%
if ~plot_approx
    ylim([60 75])
    xlim([0.01 0.1])
    legend(strcat("$d=$",num2str(d(1)),", $k_2=$ ",num2str(k2(1))),...
        strcat("$d=$",num2str(d(1)),", $k_2=$ ",num2str(k2(2))), ...
        strcat("$d=$",num2str(d(2)),", $k_2=$ ",num2str(k2(1))), ...
        strcat("$d=$",num2str(d(2)),", $k_2=$ ",num2str(k2(2))),...
        'Interpreter','latex')
end
%%
% Plot approximation.
if plot_approx
    ylim([60 75])
    xlim([0.01 0.1])
    approx_number = N_approx(L,N,L,k2(2),b/dt);
    plot(b,approx_number,'LineWidth',2)
    approx_number = N_approx(L,N,L,k2(1),b/dt);
    plot(b,approx_number,'LineWidth',2)
    legend(strcat("$d=$",num2str(d(1)),", $k_2=$",num2str(k2(1))),...
        strcat("$d=$",num2str(d(1)),", $k_2=$",num2str(k2(2))), ...
        strcat("$d=$",num2str(d(2)),", $k_2=$",num2str(k2(1))), ...
        strcat("$d=$",num2str(d(2)),", $k_2=$",num2str(k2(2))),...
        "Theoretical Fit k $= 4$",...
        "Theoretical Fit k $= 10$",...
        'Interpreter','latex')
end
%%
% Number of Cells at Steady State (Single-Cell) - function of k
'Start'

% Set up parameters for simulation.
N = 50; L = 10;
k2 = 0.5:0.5:20;
n = 100000; nb = 10000;
d = [0.00005, 0.0002]; b = [0.03, 0.07]; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 20;
ss_number = zeros(length(k2),1);
ss_number_test = zeros(tests,1);

mc=0;

f4 = figure('Visible','on');
for b_i = 1:length(b)
    for d_i = 1:length(d)
        % Run test on number
        for i = 1:length(k2)
            for j = 1:tests
        
                % If gradual death
                if ~id
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function(k2(i),k2(i),L,N,n,nb,b(b_i),d(d_i),dt,mc);
        
                % If instant death   
                else
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function_id(k2(i),k2(i),L,N,n,nb,b(b_i),d(d_i),dt,mc);
                end

                % Record data from test j
                ss_number_test(j,1) = number_ss_1;
                
                perc = 100*((i-1)*tests + j)/(length(k2)*tests);
                disp(strcat(num2str(round(perc,1)),"% finished"));
            end
        
            % Average data over all tests.
            ss_number(i) = mean(ss_number_test,'all');
        
        end
        figure(f4);
        hold on
        plot(k2,ss_number,'o-')
        xlabel('k')
        ylabel('Number of Cells')
    end
end
legend(strcat("$b=$",num2str(b(1)),", $d=$",num2str(d(1))),strcat("$b=$",...
    num2str(b(1)),", $d=$",num2str(d(2))),strcat("$b=$",num2str(b(2)),...
    ", $d=$",num2str(d(1))),strcat("$b=$",num2str(b(2)),", $d=$",num2str(d(2))),...
    'Interpreter','latex')

if plot_approx
    % Plot approximation.
    approx_number = N_approx(L,N,L,k2,b(2)/dt);
    plot(k2,approx_number)
    approx_number = N_approx(L,N,L,k2,b(1)/dt);
    plot(k2,approx_number)
    legend(strcat("$b=$",num2str(b(1)),", $d=$",num2str(d(1))), ...
        strcat("$b=$",num2str(b(1)),", $d=$",num2str(d(2))),...
        strcat("$b=$",num2str(b(2)),", $d=$",num2str(d(1))),...
        strcat("$b=$",num2str(b(2)),", $d=$",num2str(d(2))),...
        "Theoretical Fit b $= 0.07$",...
        "Theoretical Fit b $= 0.03$",...
        'Interpreter','latex')
end


% Number of Cells at Steady State (Single-Cell) - function of d.
'Start'

% Set up parameters for simulation.
N = 50; L = 10;
k2 = [4 10];
n = 100000; nb = 10000;
d = 0.00001:0.00001:0.0002; 
b = [0.03, 0.07]; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 20;
ss_number = zeros(length(d),1);
ss_number_test = zeros(tests,1);

mc=0;

f5 = figure('Visible','on');
colors = ['r','b','g','k'];
c=0;
for b_i = 1:length(b)
    % Run test on invasion speed
    for k_i = 1:length(k2)        
        c=c+1;
        for i = 1:length(d)
            for j = 1:tests
        
                % If gradual death
                if ~id
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function(k2(k_i),k2(k_i),L,N,n,nb,b(b_i),d(i),dt,mc);
        
                % If instant death   
                else
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function_id(k2(k_i),k2(k_i),L,N,n,nb,b(b_i),d(i),dt,mc);
                end
                
                % Record data from test j
                ss_number_test(j,1) = number_ss_1;
                
                perc = 100*((i-1)*tests + j)/(length(d)*tests);
                disp(strcat(num2str(round(perc,1)),"% finished"));
            end
    
            % Average data over all tests.
            ss_number(i) = mean(ss_number_test,'all');
            
        end

        figure(f5);
        hold on
        plot(d,ss_number,'Color',colors(c))
        if plot_approx
            d_N_approx = zeros(1,length(d)) + N_approx(L,N,L,k2(k_i),b(b_i)/dt);
            plot(d,d_N_approx,'--','Color',colors(c))
        end
        xlabel('d')
        ylabel('Number of Cells')
    end
end

legend(strcat("$b=$",num2str(b(1)),", $k_2=$",num2str(k2(1))),'' ...
    ,strcat("$b=$",num2str(b(1)),", $k_2=$",num2str(k2(2))),'', ...
    strcat("$b=$",num2str(b(2)),", $k_2=$",num2str(k2(1))),'', ...
    strcat("$b=$",num2str(b(2)),", $k_2=$",num2str(k2(2))),'',...
    'Interpreter','latex')

% Number of Cells at Steady State (Single-Cell) - function of b.
'Start'

% Set up parameters for simulation.
N = 50; L = 10;
k2 = [4 10];
n = 100000; nb = 10000;
d = [0.00005, 0.0002]; 
b = 0.01:0.01:0.1; dt = 1e-3;
id = 1;

% Set up parameters for testing.

tests = 20;
ss_number = zeros(length(b),1);
ss_number_test = zeros(tests,1);

mc=0;

f6 = figure('Visible','on');
for d_i = 1:length(d)
    % Run test on invasion speed
    for k_i = 1:length(k2)        
        for i = 1:length(b)
            for j = 1:tests
        
                % If gradual death
                if ~id
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function(k2(k_i),k2(k_i),L,N,n,nb,b(i),d(d_i),dt,mc);
        
                % If instant death   
                else
                    [final_time1,final_time2,winner,number_ss_1,number_ss_2,...
                        death_count_ss,variance_length,average_density,...
                        instant_speed,border_v,position_array,number_cells] = ...
                        current_invasion_function_id(k2(k_i),k2(k_i),L,N,n,nb,b(i),d(d_i),dt,mc);
                end
                
                % Record data from test j
                ss_number_test(j,1) = number_ss_1;
                
                perc = 100*((i-1)*tests + j)/(length(b)*tests);
                disp(strcat(num2str(round(perc,1)),"% finished"));
            end
    
            % Average data over all tests.
            ss_number(i) = mean(ss_number_test,'all');
            
        end

        figure(f6);
        hold on
        plot(b,ss_number,'*-')
        xlabel('b')
        ylabel('Number of Cells')
    end
end

if ~plot_approx
    legend(strcat("$b=$",num2str(b(1)),", $k_2=$ ",num2str(k2(1))),...
        strcat("$b=$",num2str(b(1)),", $k_2=$ ",num2str(k2(2))), ...
        strcat("$b=$",num2str(b(2)),", $k_2=$ ",num2str(k2(1))), ...
        strcat("$b=$",num2str(b(2)),", $k_2=$ ",num2str(k2(2))),...
        'Interpreter','latex')
end

% Plot approximation.
if plot_approx
    approx_number = N_approx(L,N,L,k2(2),b/dt);
    plot(b,approx_number)
    approx_number = N_approx(L,N,L,k2(1),b/dt);
    plot(b,approx_number)
    legend(strcat("$b=$",num2str(b(1)),", $k_2=$",num2str(k2(1))),...
        strcat("$b=$",num2str(b(1)),", $k_2=$",num2str(k2(2))), ...
        strcat("$b=$",num2str(b(2)),", $k_2=$",num2str(k2(1))), ...
        strcat("$b=$",num2str(b(2)),", $k_2=$",num2str(k2(2))),...
        "Theoretical Fit k $= 4$",...
        "Theoretical Fit b $= 10$",...
        'Interpreter','latex')
end