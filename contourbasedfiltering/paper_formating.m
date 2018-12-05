% figure
openfig('TrueContourTraj_2.fig')
h=gca
reset(h) 
axis([-6,6,-25,25])
grid on
axis square
% axis equal

plot_prop_paper_2D
saveas(gcf,'testplot','png')
pause(0.5)
close all

%%
foldername='/home/venkat/GoogleDrive/Nagavenkat_Adurthi_DRIVE/2018 papers/MOC: AAS,ACC,Journal/figures/duffprop_journal/traj+prior_surf';
% foldername='/home/venkat/GoogleDrive/Nagavenkat_Adurthi_DRIVE/2018 papers/MOC: AAS,ACC,Journal/figures/duffprop_journal/traj+prior_contour';
% foldername='/home/venkat/GoogleDrive/Nagavenkat_Adurthi_DRIVE/2018 papers/MOC: AAS,ACC,Journal/figures/duffprop_journal/prior';
% foldername='/home/venkat/GoogleDrive/Nagavenkat_Adurthi_DRIVE/2018 papers/MOC: AAS,ACC,Journal/figures/duffprop_journal/post';


clc
files=dir(foldername);

parfor i=1:length(files)
    if files(i).isdir==0
        if contains(files(i).name,'.fig')==1
            fpath = [files(i).folder,'/',files(i).name];
            customplotmodifierfunc(fpath);
        end
    end
end