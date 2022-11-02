function xvec_new = border_test(x)
% Function to check if cell boundary positions overlap i.e overshoot. Fix
% if so.
N = length(x)-1;
x_old = x;
for i=2:N
    if x_old(i) >= x_old(i+1)
        x_old(i) = x_old(i+1)*0.999;
    end
    if x_old(i) <= x_old(i-1)
        x_old(i) = x_old(i-1)*1.001;
    end
end
xvec_new = x_old;