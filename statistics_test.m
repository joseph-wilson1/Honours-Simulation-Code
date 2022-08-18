%% Invasion Speed
'Start'
clear all;
k1 = 5; N = 50; L = 10;
k2 = 0.5:0.5:5;
%n = 10000; nb = 1000;
n = 1000000; nb = 10000;
d = 0.002; b = 0.7; dt = 1e-2;
tests = 10;
speeds = zeros(length(k2),tests);
death_count = zeros(length(k2),tests);
speed_ave = zeros(length(k2),1);
death_ave = zeros(length(k2),1);
for i = 1:length(k2)
    for j = 1:tests
        [final_time1,final_time2,density,winner,~,SS_death,~] = current_invasion_function(k1,k2(i),L,N,n,nb,b,d,dt,1);
        %[final_time1,final_time2,density,winner,~,~] = current_invasion_function_id(k1,k2(i),L,N,n,nb,b,d,dt,1);
        if winner
            if final_time1==0
                speeds(i,j) = L/2/final_time2;
            elseif final_time2==0
                speeds(i,j) = -L/2/final_time1;
            end
        else
            speeds(i,j) = 0;
        end
        death_count(i,j) = SS_death;
    end
    speed_ave(i) = mean(speeds(i,:));
    death_ave(i) = mean(death_count(i,:));
end

%%
subplot(2,1,1)
plot(k2,speed_ave,'o--')
title("Invasion Speed vs k_2  -> k_1 = 5")
xlabel("k_2")
ylabel("Invasion Speed")
subplot(2,1,2)
plot(k2,death_ave,'o--')
title("Death Count at Steady State")
xlabel("k_2")
ylabel('Death Count')


%% Steady State Density
clear all;
N = 100; L = 10;
k2 = 5;
n = 10000; nb = 5000;
d = 0.002; b = 0.7; dt = 1e-2;
k = 0.5:0.5:10;
tests = 5;
density_vec = zeros(1,tests);
SS_density_ave = zeros(length(k),1);
variance_vec = zeros(1,tests);
variance_ave = zeros(length(k),1);
for i = 1:length(k)
    for j = 1:tests
        [~,~,~,~,SS_density,~,variance_length] = current_invasion_function(k(i),k2,L,N,n,nb,b,d,dt,0);
        density_vec(j) = SS_density;
        variance_vec(j) = variance_length;
    end
    variance_ave(i) = mean(variance_vec,'all');
    SS_density_ave(i) = mean(density_vec,'all');
end

%%
plot(k,variance_ave,'o--')
title(strcat("Average Variance of Length (Steady State) vs k for N = ",num2str(N),"cells"))
xlabel('k')
ylabel('Variance of Length')
legend('d = 0.002, b = 0.7')

%% Steady State Density (function of k)
clear all;
N = 50; L = 10;
k1 = 0.5:0.5:20; 
k2 = 10;
n = 100000; nb = 50000;
d = [0.00005,0.001]; dt = 1e-3; % d = [0.05,0.5]
b = [0.03 0.07]; % b = [30,70]
tests = 5;
density_vec = zeros(1,tests);
% death_vec = zeros(1,tests);
SS_density_ave = zeros(length(k1),1);
% SS_death_ave = zeros(length(k1),1);
for l=1:length(d)
    l
    for k = 1:length(b)
        for i = 1:length(k1)
            for j = 1:tests
                [~,~,~,~,SS_density,death_count,~] = current_invasion_function(k1(i),k2,L,N,n,nb,b(k),d(l),dt,0);
                density_vec(j) = SS_density;
%                 death_vec(j) = mean(death_count(1,nb:end));
            end
%             SS_death_ave(i) = mean(death_vec,'all');
            SS_density_ave(i) = mean(density_vec,'all');
        end
        figure(1)
        hold on
        plot(k1,SS_density_ave,'o-')
        xlabel('k')
        ylabel('Number of Cells')
        ylim([50 75])
        xlim([k1(1) k1(end)])
    end
end
legend("b=0.03,d=0.00005","b=0.07,d=0.00005","b=0.03,d=0.001","b=0.07,d=0.001")

%%
k1 = 0.5:0.5:20; 
f2 = @(x,b) 50/12*(14+exp(-3*x/b)+3*exp(-x/b));
y1 = f2(k1,70);
y2 = f2(k1,30);
hold on
plot(k1,y1)
plot(k1,y2)
legend("b=0.03,d=0.00005","b=0.07,d=0.00005","b=0.03,d=0.001","b=0.07,d=0.001",...
"Theoretical Fit b = 0.07","Theoretical Fit b = 0.03")

%% Steady State Density (function of b)
clear all;
N = 50; L = 10;
k1 = [4,10]; k2 = 5;
n = 100000; nb = 50000;
d = [0.00005,0.001]; dt = 1e-3;
b = 0.01:0.005:0.1;
tests = 3;
density_vec = zeros(1,tests);
death_vec = zeros(1,tests);
SS_density_ave = zeros(length(b),1);
SS_death_ave = zeros(length(b),1);
for l=1:length(d)
    l
    for k = 1:length(k1)
        k
        for i = 1:length(b)
            for j = 1:tests
                [~,~,~,~,SS_density,~,~] = current_invasion_function(k1(k),k2,L,N,n,nb,b(i),d(l),dt,0);
                density_vec(j) = SS_density;
%                 death_vec(j) = mean(death_count(1,nb:end));
            end
%             SS_death_ave(i) = mean(death_vec,'all');
            SS_density_ave(i) = mean(density_vec,'all');
        end
        figure(1)
        hold on
        plot(b,SS_density_ave,'*-')
        xlabel('b')
        ylabel('Number of Cells')
        ylim([50 75])
        xlim([b(1) b(end)])
    end
end
legend("k=4,d=0.00005","k=10,d=0.00005","k=4,d=0.001","k=10,d=0.001")

%%
b = 5:5:100;
f1 = @(x,k) 50/12*(14+exp(-3*k./x)+3*exp(-k./x));
y1 = f1(b,4);
y2 = f1(b,10);
hold on
plot(b/1000,y1,'LineWidth',2)
plot(b/1000,y2,'LineWidth',2)
legend("k=4,d=0.00005","k=10,d=0.00005","k=4,d=0.001","k=10,d=0.001",...
"Theoretical Fit k = 4","Theoretical Fit k = 10")

%% Steady State Density (function of d)
clear all;
N = 50; L = 10;
k1 = [4,10]; k2 = 5;
n = 100000; nb = 50000;
d = 0.00005:0.00005:0.001; dt = 1e-3;
b = [0.03 0.07];
tests = 10;
density_vec = zeros(1,tests);
death_vec = zeros(1,tests);
SS_density_ave = zeros(length(d),1);
SS_death_ave = zeros(length(d),1);
for l=1:length(b)
    for k = 1:length(k1)
        for i = 1:length(d)
            for j = 1:tests
                [~,~,~,~,SS_density,death_count,~] = current_invasion_function(k1(k),k2,L,N,n,nb,b(l),d(i),dt,0);
                density_vec(j) = SS_density;
                death_vec(j) = mean(death_count(1,nb:end));
            end
            SS_death_ave(i) = mean(death_vec,'all');
            SS_density_ave(i) = mean(density_vec,'all');
        end
        figure(1)
        hold on
        plot(d,SS_density_ave,'+-')
        xlabel('d')
        ylabel('Number of Cells')
        ylim([50 75])
        xlim([d(1) d(end)])
    end
end
legend("k=4,b=0.03","k=10,b=0.03","k=4,b=0.07","k=10,b=0.07")

%% Dead Cells
clear all;
N = 100; L = 10;
k2 = 5;
n = 10000; nb = 5000;
d = 0.002; b = 0.7; dt = 1e-2;
k = 0.5:0.5:10;
for i = 1:length(k)
    [~,~,~,~,~,death_count] = current_invasion_function(k(i),k2,L,N,n,nb,b,d,dt,0);
    hold on
    subplot(4,5,i)
    histogram(death_count,'Normalization','probability')
    title(strcat('K = ',num2str(k(i))))
    xlabel('Dead Cell Count')
    ylabel('Probability')
end
sgtitle(strcat('Normalized histogram of Dead Cell Count before SS vs k - with N =',num2str(N),'cells'))

