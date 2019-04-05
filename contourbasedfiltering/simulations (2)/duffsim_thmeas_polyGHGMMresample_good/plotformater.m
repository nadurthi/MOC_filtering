foldername ='simulations/duffsim_thmeas_polyGHGMMresample_good/post';

clc
files=dir(foldername);

for i=1:length(files)
    if files(i).isdir==0
        if contains(files(i).name,'.fig')==1 && contains(files(i).name,'PZ')==1
            fpath = [files(i).folder,'/',files(i).name];
            openfig(fpath)
             if contains(files(i).name,'Surf')==1
%                  view([-17,65])
             end
             grid
%             customplotmodifierfunc(fpath);
            plot_prop_paper_2D
            axis square
            saveas(gcf,fpath(1:end-4),'png')
            % saveas(gcf,fpath(1:end-4),'fig')
            pause(0.1)
            close
        end
    end
end
