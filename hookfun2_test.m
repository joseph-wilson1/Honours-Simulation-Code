'Start'
L = 10; x0 = 0; N = 10;
k0 = 10;
k = zeros(1,N)+k0;
a0 = L/N; 
a = zeros(1,N)+a0;
dt = 1e-2; eta = 1;
m = dt/eta;
xvec = linspace(x0,L,N+1);
xvec
xvec_new = hookfun2(xvec,k,a,dt,eta)

x = xvec;
x9 = m*k(9)*x(10) + (1-m*k(9) - m*k(8))*x(9) + m*k(8)*x(8) - m*k(9)*a(9)...
    +m*k(8)*a(8);
x10 = m*k(10)*x(11) + (1-m*k(10) - m*k(9))*x(10) + m*k(9)*x(9) - m*k(10)*a(10)...
    +m*k(9)*a(9);

x3 = m*k(3)*xvec(3)+(1-m*(k(3)+k(3-1)))*xvec(3-1)+m*k(3-1)*xvec(3-2)...
    -m*(k(3)*a(3)-k(3-1)*a(3-1));