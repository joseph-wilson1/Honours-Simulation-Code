b0 = 1;
n = 100;
dt = 1/n;
b = b0*dt;
x = 0:0.1:20;


y = dt*geopdf(x,b);
y2 = exppdf(x,b0) % y2 = y*dt

hold off
plot(x,y)
hold on
plot(x,y2)
legend("geo",'exp')