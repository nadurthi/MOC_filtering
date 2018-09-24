%% marginal 0I 2D plots

plotmargs=1;

if plotmargs == 1
    if isempty(Xmctest)==0
        ind = sqrt(sum(Xnmctest.^2,2))<2.3;
        Xnmctest=Xnmctest(ind,:);
    end
    
    [Xx,Xy]=meshgrid(linspace(-2,2,50),linspace(-2,2,50) );
    % Xp=[reshape(Xx,625,1),reshape(Xy,625,1)];
    margprobs = zeros(size(Xx));
    %     margprobs_cell = cell(size(Xx,1),1);
    for i=1:size(Xx,1)
        for j=1:size(Xx,2)
            margprobs(i,j) = pdfnorm.func([Xx(i,j),Xy(i,j)]);
        end
    end
    % get the true points and their probs
    margprobs_true = zeros(size(Xx));
    Xx_true = zeros(size(Xx));
    Xy_true = zeros(size(Xx));
    for i=1:size(Xx,1)
        for j=1:size(Xx,2)
            g=pdfnorm.normX2trueX([Xx(i,j),Xy(i,j)]);
            Xx_true(i,j) = g(1);
            Xy_true(i,j) = g(2);
            margprobs_true(i,j) = pdfnorm.normprob2trueprob( margprobs(i,j) );
        end
    end
    
    
    figure(1)
    contour(Xx,Xy,margprobs,15)
    grid on
    box off
    hold on
    if isempty(Xmctest)==0
        plot(Xnmctest(:,1),Xnmctest(:,2),'r.')
    end
    if isempty(Xtruth)==0
        if plotsconf.fig1.plottruth == true
            plot(Xntruth(:,1),Xntruth(:,2),'k*','linewidth',2)
        end
    end
    
    %     plot(Xt(:,1),Xt(:,2),'g*')
    title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
    xlabel('x_1')
    ylabel('x_2')
    axis equal
    axis square
    hold off
    saveas(gcf,[plotsconf.plotfolder,'/NormContour_',plotsconf.nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[plotsconf.plotfolder,'/NormContour_',plotsconf.nametag,'_',num2str(Tk)],'fig')
    
    figure(2)
    surf(Xx,Xy,margprobs,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
    camlight right; lighting phong
    alpha 0.4
    hold on
    if isempty(Xmctest)==0
        plot(Xnmctest(:,1),Xnmctest(:,2),'r.')
    end
    if isempty(Xtruth)==0
        if plotsconf.fig2.plottruth == true
            plot(Xntruth(:,1),Xntruth(:,2),'k*','linewidth',2)
        end
    end
    
    %     plot(Xt(:,1),Xt(:,2),'g*')
    title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
    xlabel('x_1')
    ylabel('x_2')
    axis equal
    axis square
    hold off
    saveas(gcf,[plotsconf.plotfolder,'/NormSurf_',plotsconf.nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[plotsconf.plotfolder,'/NormSurf_',plotsconf.nametag,'_',num2str(Tk)],'fig')
    
    figure(3)
    if plotsconf.fig3.holdon == true
        hold on
    end
    [~,h]=contour(Xx_true,Xy_true,margprobs_true ,15);
    grid on
    box off
    h.ContourZLevel =plotsconf.fig3.contourZshift;
    view([26,43])
    hold on
    if isempty(Xmctest)==0
        plot3(Xmctest(:,1),Xmctest(:,2),repmat(plotsconf.fig3.contourZshift,size(Xmctest,1),1),'r.')
    end
    if isempty(Xtruth)==0
        if plotsconf.fig3.plottruth == true
            plot3(Xtruth(:,1),Xtruth(:,2),plotsconf.fig3.contourZshift,'k*','linewidth',2)
        end
    end
    if isempty(plotsconf.fig3.plotmeas)==0
        zz = plotsconf.fig3.plotmeas;
        plot3(zz(1),zz(2),plotsconf.fig3.contourZshift,'bs','linewidth',2)
    end
    %     plot(Xt(:,1),Xt(:,2),'g*')
    title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
    xlabel('x_1')
    ylabel('x_2')
    axis equal
    axis square
    hold off
    saveas(gcf,[plotsconf.plotfolder,'/TrueContour_',plotsconf.nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[plotsconf.plotfolder,'/TrueContour_',plotsconf.nametag,'_',num2str(Tk)],'fig')
    
    figure(4)
    if plotsconf.fig4.holdon == true
        hold on
    end
    surf(Xx_true,Xy_true,margprobs_true,'FaceColor',plotsconf.fig4.surfcol,'EdgeColor','none','FaceAlpha',0.7);
    camlight right; lighting phong
    alpha 0.7
    hold on
    if isempty(Xmctest)==0
        plot(Xmctest(:,1),Xmctest(:,2),'r.')
    end
    if isempty(Xtruth)==0
        if plotsconf.fig4.plottruth == true
            plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
        end
    end
    if isempty(plotsconf.fig4.plotmeas)==0
        zz = plotsconf.fig4.plotmeas;
        plot(zz(1),zz(2),'bs','linewidth',2)
    end
    
    %     plot(Xt(:,1),Xt(:,2),'g*')
    title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
    xlabel('x_1')
    ylabel('x_2')
    axis equal
    axis square
    hold off
    saveas(gcf,[plotsconf.plotfolder,'/TrueSurf_',plotsconf.nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[plotsconf.plotfolder,'/TrueSurf_',plotsconf.nametag,'_',num2str(Tk)],'fig')
    
end