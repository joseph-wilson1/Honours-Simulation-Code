function plot_density(density,n,dt,multiple_cell,k1,k2)
% FUNCTION NOT IN USE.
    figure(3)
    clf;
    if multiple_cell
        plot(0:dt:dt*n,density(1,:))
        hold on
        plot(0:dt:dt*n,density(2,:))
        title("Number of cells vs Time")
        xlabel("time")
        ylabel("number of cells")
        legend(strcat("Cell 1, k = ",num2str(k1)),strcat("Cell 2, k = ",num2str(k2)))
    else
        plot(0:dt:dt*n,density(1,:))
        title("Number of cells vs Time")
        xlabel("time")
        ylabel("number of cells")
    end
end