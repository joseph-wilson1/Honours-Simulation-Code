L = 10; x0 = 0; N = 10;
k0 = 10; 
border = N/2+1;
k = zeros(1,N)+k0;
dead = zeros(1,N);
k(1) = 9;
dead(6) = 1;
dead(1)=1;
a0 = L/N; length_tol = a0/10; 
a = zeros(1,N)+a0;
a(5) = 2;
dt = 1e-2; eta = 1;
xvec = linspace(x0,L,N+1);
xvec(2) = xvec(1) + (a0/11);
xvec(6) = xvec(7) - (a0/11);
xvec
k
dead
a
[xvec_new,k_new,a_new,dead_new,~,~] = length_test(xvec,k,a,dead,length_tol,border)

