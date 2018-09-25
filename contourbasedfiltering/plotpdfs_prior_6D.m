function plotpdfs_prior_6D(states,Tk,pdfnorm,X,probs,mquadf,Pquadf,Xmc,Xtruth,saveprops)
nametag='prior';

%%
dim=size(X,2);
states2keep = states;
states2remove = 1:dim;
states2remove(states2keep)=[];


Xn = pdfnorm.transForms.trueX2normX(X);
pn = pdfnorm.transForms.trueprob2normprob(probs);
GMMfitter = GMMFitDataSet(Xn,pn);
GMM = GMMfitter.FitGMM_kmeans_optimwt(3);



if isempty(Xmc)==0
    Xmcnorm = pdfnorm.transForms.trueX2normX(Xmc) ;
end
if isempty(Xtruth)==0
    Xtruthnorm = pdfnorm.transForms.trueX2normX(Xtruth) ;
end
if isempty(mquadf)==0
    [x,w] = UT_sigmapoints(mquadf,Pquadf,2);
    x = pdfnorm.transForms.trueX2normX(x) ;
    [mquadfnorm,Pquadfnorm]=MeanCov(x,w);
    mquadfnorm=mquadfnorm(:);
end

[mX,PX] = MeanCov(Xn,pn/sum(pn));

[Xx,Xy]=meshgrid(linspace(-2,2,25),linspace(-2,2,25) );


pdfprobs_norm_cell = cell(size(Xx,1),1);
QuadFilprobs_norm_cell = cell(size(Xx,1),1);
parfor i=1:size(Xx,1)
    pdfprobs_norm = zeros(size(Xx,2),1);
    QuadFilprobs_norm = zeros(size(Xx,2),1);
    Xpoint = zeros(1,dim);
    for j=1:size(Xx,2)
        Xpoint(states2keep) = [Xx(i,j),Xy(i,j)];
        Xpoint(states2remove) = mX(states2remove);
%         keyboard
%         pdfprobs_norm(j) = pdfnorm.func(Xpoint);
%         QuadFilprobs_norm(j) = mvnpdf(Xpoint,mquadfnorm',Pquadfnorm);
        
        pdfprobs_norm(j) = marginalize_exp_pdf_modf([Xx(i,j),Xy(i,j)],states2remove,pdfnorm,X,probs,mX,PX,GMM,'GMM_MC');
        QuadFilprobs_norm(j) = mvnpdf([Xx(i,j),Xy(i,j)],mquadfnorm(states)',Pquadfnorm(states,states));
    end
    pdfprobs_norm_cell{i}=pdfprobs_norm;
    QuadFilprobs_norm_cell{i}=QuadFilprobs_norm;
end

pdfprobs_norm = zeros(size(Xx));
QuadFilprobs_norm = zeros(size(Xx));
for i=1:size(Xx,1)
    for j=1:size(Xx,2)
        pdfprobs_norm(i,j) = pdfprobs_norm_cell{i}(j);
        QuadFilprobs_norm(i,j) = QuadFilprobs_norm_cell{i}(j);
    end
end

% % get the true points and their probs
% [Xx_true,Xy_true]=meshgrid(linspace(-2,2,50),linspace(-2,2,50) );
% pdfprobs_true = zeros(size(Xx));
% QuadFilprobs_true = zeros(size(Xx));
% Xx_true = zeros(size(Xx));
% Xy_true = zeros(size(Xx));
% for i=1:size(Xx,1)
%     for j=1:size(Xx,2)
%         g=pdfnorm.transForms.normX2trueX([Xx(i,j),Xy(i,j)]);
%         Xx_true(i,j) = g(1);
%         Xy_true(i,j) = g(2);
%         pdfprobs_true(i,j) = pdfnorm.transForms.normprob2trueprob( pdfprobs_norm(i,j) );
%         QuadFilprobs_true(i,j) = mvnpdf([Xx_true(i,j),Xy_true(i,j)],mquadf(:)',Pquadf);
%     end
% end

%%
figure(1)
contour(Xx,Xy,pdfprobs_norm,15)
grid on
box off
hold on
if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(2)
surf(Xx,Xy,pdfprobs_norm,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
% plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruthnorm)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormSurf_',nametag,'_',num2str(Tk)],'fig')
end
%%
% figure(3)
% [~,h]=contour(Xx_true,Xy_true,pdfprobs_true ,15);
% grid on
% box off
% % h.ContourZLevel = 1;
% % view([26,43])
% hold on
% if isempty(Xmc)==0
%     plot(Xmc(:,1),Xmc(:,2),'r.')
% end
% if isempty(Xtruth)==0
%     plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
% end
% 
% % title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
% xlabel('x_1')
% ylabel('x_2')
% axis equal
% axis square
% hold off
% if saveprops.saveit==1
%     saveas(gcf,[saveprops.plotfolder,'/TrueContour_',nametag,'_',num2str(Tk)],'png')
%     saveas(gcf,[saveprops.plotfolder,'/TrueContour_',nametag,'_',num2str(Tk)],'fig')
% end
% %%
% figure(4)
% surf(Xx_true,Xy_true,pdfprobs_true,'FaceColor','r','EdgeColor','none','FaceAlpha',0.7);
% camlight right; lighting phong
% alpha 0.7
% hold on
% plot3(X(:,1),X(:,2),probs,'bo')
% 
% if isempty(Xmc)==0
%     plot(Xmc(:,1),Xmc(:,2),'r.')
% end
% if isempty(Xtruth)==0
%     plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
% end
% 
% % title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
% xlabel('x_1')
% ylabel('x_2')
% axis equal
% axis square
% hold off
% if saveprops.saveit==1
%     saveas(gcf,[saveprops.plotfolder,'/TrueSurf_',nametag,'_',num2str(Tk)],'png')
%     saveas(gcf,[saveprops.plotfolder,'/TrueSurf_',nametag,'_',num2str(Tk)],'fig')
% end

%%
% -------------------- Quad Filter Plots------------------------------------------------------------------------
%%
%%
figure(5)
contour(Xx,Xy,QuadFilprobs_norm,15)
grid on
box off
hold on
if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(6)
surf(Xx,Xy,QuadFilprobs_norm,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
% plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormSurf_',nametag,'_',num2str(Tk)],'fig')
end
%%
% figure(7)
% [~,h]=contour(Xx_true,Xy_true,QuadFilprobs_true ,15);
% grid on
% box off
% % h.ContourZLevel = 1;
% % view([26,43])
% hold on
% if isempty(Xmc)==0
%     plot(Xmc(:,1),Xmc(:,2),'r.')
% end
% if isempty(Xtruth)==0
%     plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
% end
% 
% % title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
% xlabel('x_1')
% ylabel('x_2')
% axis equal
% axis square
% hold off
% if saveprops.saveit==1
%     saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueContour_',nametag,'_',num2str(Tk)],'png')
%     saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueContour_',nametag,'_',num2str(Tk)],'fig')
% end
% %%
% figure(8)
% surf(Xx_true,Xy_true,QuadFilprobs_true,'FaceColor','r','EdgeColor','none','FaceAlpha',0.7);
% camlight right; lighting phong
% alpha 0.7
% hold on
% plot3(X(:,1),X(:,2),probs,'bo')
% 
% if isempty(Xmc)==0
%     plot(Xmc(:,1),Xmc(:,2),'r.')
% end
% if isempty(Xtruth)==0
%     plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
% end
% 
% % title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
% xlabel('x_1')
% ylabel('x_2')
% axis equal
% axis square
% hold off
% if saveprops.saveit==1
%     saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueSurf_',nametag,'_',num2str(Tk)],'png')
%     saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueSurf_',nametag,'_',num2str(Tk)],'fig')
% end

end