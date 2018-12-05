close all
ff = 'simulations/duffsim_rthmeas_Rr1th2deg_polyGHSampler';
mkdir([ff,'/traj+priorpost_contour'])
% figure(1)
% plot(Xtruplot(:,1),Xtruplot(:,2),'k--')
% hold on
for k=2:time.Ntsteps
    Xmctest = zeros(size(XMC,1),model.fn);
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
    Xmctest = [];
    
    truecontpriorfig= [ff,'/prior/TrueContour_prior_',num2str(k),'.fig'];
    truecontpostfig= [ff,'/post/TrueContour_post_',num2str(k),'.fig'];
    openfig(truecontpriorfig)
    h = gcf;
    axesObjs = get(h, 'Children');
    dataObjsprior = get(axesObjs, 'Children');
    
    pause(1)
    
    openfig(truecontpostfig)
    h = gcf;
    axesObjs = get(h, 'Children');
    dataObjspost = get(axesObjs, 'Children');
    
    figure(301)
    contour(dataObjsprior(2).XData,dataObjsprior(2).YData,dataObjsprior(2).ZData,dataObjsprior(2).LevelList,'r')
    hold on
    contour(dataObjspost(2).XData,dataObjspost(2).YData,dataObjspost(2).ZData,dataObjspost(2).LevelList,'b')
    plot(Xtruplot(:,1),Xtruplot(:,2),'k--')
%     plot(Xmctest(:,1),Xmctest(:,2),'ro')
    axis([-10,10,-40,40])
    title(['k = ',num2str(k)])
%     plot_prop_paper
    hold off
    
    saveas(gcf,[ff,'/traj+priorpost_contour','/TrueContourTrajNOMC_',num2str(k)],'png')
    saveas(gcf,[ff,'/traj+priorpost_contour','/TrueContourTrajNOMC_',num2str(k)],'fig')
    
    pause(1)
    

    figure(302)
    surf(dataObjsprior(2).XData,dataObjsprior(2).YData,dataObjsprior(2).ZData,'FaceColor','b','EdgeColor','none','FaceAlpha',0.7)
    camlight right; lighting phong
    alpha 0.7
    box off
    grid on

    hold on
    
    surf(dataObjspost(2).XData,dataObjspost(2).YData,dataObjspost(2).ZData,'FaceColor','r','EdgeColor','none','FaceAlpha',0.7)
    plot(Xtruplot(:,1),Xtruplot(:,2),'k--')
    camlight right; lighting phong
    alpha 0.7
    hold off
    
    %     plot(Xmctest(:,1),Xmctest(:,2),'ro')
    axis([-10,10,-40,40])
    title(['k = ',num2str(k)])
    view([0,56])
%     plot_prop_paper
    
    
    
    saveas(gcf,[ff,'/traj+priorpost_contour','/TrueSurfTrajNOMC_',num2str(k)],'png')
    saveas(gcf,[ff,'/traj+priorpost_contour','/TrueSurfTrajNOMC_',num2str(k)],'fig')

    pause(1)
end

