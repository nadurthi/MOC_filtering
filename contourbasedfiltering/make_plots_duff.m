close all
ff = 'duffsim_linearmeas/prior';
mkdir('duffsim_linearmeas/traj+prior_contour')
% figure(1)
% plot(Xtruplot(:,1),Xtruplot(:,2),'k--')
% hold on
for k=2:51
    Xmctest = zeros(size(XMC,1),model.fn);
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
    
    truesurffig= [ff,'/TrueContour_prior_',num2str(k),'.fig'];
    openfig(truesurffig)
    hold on
    plot(Xtruplot(:,1),Xtruplot(:,2),'k--')
%     plot(Xmctest(:,1),Xmctest(:,2),'ro')
    axis([-10,10,-40,40])
    title(['k = ',num2str(k)])
    plot_prop_paper
    pause(1)
%     saveas(gcf,['duffsim_linearmeas/traj+prior_contour','/TrueContourTrajNOMC_',num2str(k)],'png')
%     saveas(gcf,['duffsim_linearmeas/traj+prior_contour','/TrueContourTrajNOMC_',num2str(k)],'fig')
%     f1 = openfig(filename1);
%     ax = gca;
%     x1 = ax.XData;
%     y1 = ax.YData;
%     figure(1)
%     plot(x1, y1, 'b-', 'LineWidth', 2);

end

