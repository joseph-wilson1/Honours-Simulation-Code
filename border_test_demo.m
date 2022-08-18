L = 10; x0 = 0; N = 50;
k0 = 10;
k = zeros(1,N)+k0;
a0 = L/N; 
a = zeros(1,N)+a0;
dt = 1e-2; eta = 1;
xvec = linspace(x0,L,N+1);
xvec(2) = x0-(a0/10);
xvec(end-1) = L + (a0/10);
xvec(25) = xvec(26)+(a0/10);
xvec(30) = xvec(29)-(a0/10);
% y = zeros(1,N+1);
% plot(xvec,y,'o--')
xvec
xvec_new = border_test(xvec)