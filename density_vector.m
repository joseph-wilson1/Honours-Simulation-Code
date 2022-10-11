function d_vec = density_vector(xvec)
%%% Function to calculate the density of all cells in a simulation, where
%%% density = 1/length.
N = length(xvec)-1;
d_vec = zeros(1,N);
for i=1:N
    d_vec(i) = 1/li(xvec,i);
end