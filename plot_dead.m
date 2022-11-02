function plot_dead(death_count,n,dt,k1,k2,multiple_cell)
%Plot Dead Cells. NOT IN USE.
figure(4)
clf;
if multiple_cell
    plot(0:dt:dt*n,death_count(1,:),'o')
    hold on
    plot(0:dt:dt*n,death_count(2,:),'o')
    title("Number of Dead Cells vs Time")
    xlabel("time")
    ylabel("number of Dead Cells")
    legend(strcat("Cell 1, k = ",num2str(k1)),strcat("Cell 2, k = ",num2str(k2)))
else
    plot(0:dt:dt*n,death_count(1,:),'o')
    title("Number of Dead Cells vs Time")
    xlabel("time")
    ylabel("number of Dead Cells")
end