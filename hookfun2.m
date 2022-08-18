function x_new = hookfun2(x,k,a,dt,eta)
%%% INPUTS: x = vector of old cell boundary positions, k = vector of spring
%%% constants, a = vector of resting spring lengths, dt = time-step, eta =
%%% viscosity constant
m = dt/eta;
k = [k 0];
a = [a 0];
x_new = m*k.*circshift(x,-1) + (1 - m*k - m*circshift(k,1)).*x ...
    + m*circshift(k,1).*circshift(x,1) - m*k.*a ...
    + m*circshift(k,1).*circshift(a,1);
x_new(1) = x(1);
x_new(end) = x(end);
end