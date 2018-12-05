function customplotmodifierfunc(fpath)
% fpath = [files(i).folder,'/',files(i).name];
openfig(fpath)
%             h=gca;
%             reset(h) 
%             axis([-6,6,-25,25])
%             grid on
%             axis square
plot_prop_paper_2D
saveas(gcf,fpath(1:end-4),'png')
% saveas(gcf,fpath(1:end-4),'fig')
pause(0.1)
close