close all
clear
ff = 'simulations/duffsim_prop_nonormalization_verygood';
load([ff,'/sim1.mat'])
mkdir([ff,'/traj+prior_contour'])
mkdir([ff,'/traj+prior_surf'])
% figure(1)
% plot(Xtruplot(:,1),Xtruplot(:,2),'k--')
% hold on
for k=2:51
    Xmctest = zeros(size(XMC,1),model.fn);
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
    
    truesurffig= [ff,'/prior/TrueContour_prior_',num2str(k),'.fig'];
    openfig(truesurffig)
    h = gcf;
    axesObjs = get(h, 'Children');
    dataObjsprior = get(axesObjs, 'Children');
    
    hold on
    plot(Xtruplot(:,1),Xtruplot(:,2),'k--')
    axis([-10,10,-40,40])
    title(['k = ',num2str(k)])
    plot_prop_paper
    
    pause(2)
    
%     truesurffig= [ff,'/prior/TrueSurf_prior_',num2str(k),'.fig'];
%     openfig(truesurffig)
%     h = gcf;
%     axesObjs = get(h, 'Children');
%     dataObjsprior = get(axesObjs, 'Children');
    
    pause(2)
    PPP = dataObjsprior{10}(2);
    
    figure(33)
    surf(PPP.XData,PPP.YData,PPP.ZData,'FaceColor','b','EdgeColor','none','FaceAlpha',0.7)
    camlight right; lighting phong
    alpha 0.7
    box off
    grid on
    
    hold on
    plot(Xtruplot(:,1),Xtruplot(:,2),'k--')
%     plot(Xmctest(:,1),Xmctest(:,2),'ro')
    axis([-10,10,-40,40])
    title(['k = ',num2str(k)])
    plot_prop_paper
    view([-16,44])
    
    pause(1)
    figure(33)
    saveas(gcf,[ff,'/traj+prior_surf','/TrueSurfTrajNOMC_',num2str(k)],'png')
    saveas(gcf,[ff,'/traj+prior_surf','/TrueSurfTrajNOMC_',num2str(k)],'fig')
    
    pause(1)
    close all
end

