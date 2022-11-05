h = gcf;
xlabel("$b$",Interpreter="latex")
ylabel("Number of Cells",Interpreter="latex")
filename = "birth_id_d_00005_0002_k_4_10_dt_001_n_10_5_approx";
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,filename,'-dpdf','-r0')
