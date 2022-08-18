%% Homogeneous Population

'Start'
clear all;
k1 = 10; k2 = 4;
N = 50; L = 10;
n = 50000; nb = 5000;
d = 0.0001; b = 0.7; dt = 1e-2;

[~,~,SS2,~,SS1,~] = current_invasion_function(k1,k2,L,N,n,nb,b,d,dt,1);

%%
times = 0:n;
times = times*dt;
plot(times,SS1)
plot(times,SS2)
xlabel('t')
ylabel('Number of cells')
xlim([0 100])