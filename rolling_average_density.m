function ave = rolling_average_density(xvec)
% Function to find rolling average density profile. (NOT IN USE).
N = length(xvec)-1;
d_vec = density_vector(xvec);
ave = zeros(1,N);
ave(1) = d_vec(1); ave(end) = d_vec(end);
for i=2:(length(d_vec)-1)
    ave(1,i) = (d_vec(i-1) + d_vec(i) + d_vec(i+1))/3;
end

