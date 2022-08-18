function plot_cell(xvec,border,k1,k2,x0,L,multiple_cell,j,dt)

    N = length(xvec)-1;
    if multiple_cell
        hold on
        Ys = zeros(length(xvec(1:border)),1) + j*dt;
        plot(xvec(1:border),Ys,'b.')
%         title(strcat("Heterogeneous Cell Population - k1 = ",num2str(k1)," k2 = ",num2str(k2)))
        xlabel("x")
        ylabel("t")
        xlim([x0 L])
        %Replot second population of cells
        Yd = zeros(length(xvec(border:N+1)),1) + j*dt;
        plot(xvec(border:N+1),Yd,'r.')
        xlim([x0 L])
        drawnow
    else
        hold on
        Ys = zeros(length(xvec),1) + j*dt;
        plot(xvec,Ys,'b.')
%         title(strcat("Homogenous Cell Population - k = ",num2str(k1)))
        xlabel("x")
        ylabel("t")
        xlim([x0 L])
        drawnow
        end
end