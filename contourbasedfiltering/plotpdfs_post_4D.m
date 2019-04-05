function plotpdfs_post_4D(states,Tk,pdfnormprior,pdfnorm,model,X,probs,mquadf,Pquadf,Xmc,Xtruth,zk,saveprops)
nametag='post';

%%
dim=size(X,2);
states2keep = states;
states2remove = 1:dim;
states2remove(states2keep)=[];


Xn = pdfnorm.transForms.trueX2normX(X);
pn = pdfnorm.transForms.trueprob2normprob(probs);
% GMMfitter = GMMFitDataSet(Xn,pn);
% GMM = GMMfitter.FitGMM_kmeans_optimwt(3);



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

[mXn,PXn] = MeanCov(Xn,pn/sum(pn));

[Xx,Xy]=meshgrid(linspace(-2,2,50),linspace(-2,2,50) );


pdfprobs_norm_cell = cell(size(Xx,1),1);
priorpdfprobs_norm_cell = cell(size(Xx,1),1);
QuadFilprobs_norm_cell = cell(size(Xx,1),1);

for i=1:size(Xx,1)
    pdfprobs_norm = zeros(size(Xx,2),1);
    priorpdfprobs_norm = zeros(size(Xx,2),1);
    QuadFilprobs_norm = zeros(size(Xx,2),1);
    Xpoint = zeros(size(Xx,2),dim);
    Xpointquad = zeros(size(Xx,2),dim);
    for j=1:size(Xx,2)
        Xpoint(j,states2keep) = [Xx(i,j),Xy(i,j)];
        Xpoint(j,states2remove) = mXn(states2remove);
        
        Xpointquad(j,states2keep) = [Xx(i,j),Xy(i,j)];
        Xpointquad(j,states2remove) = mquadfnorm(states2remove);
    end
    pdfprobs_norm = pdfnorm.func(Xpoint);
    QuadFilprobs_norm = mvnpdf(Xpointquad,mquadfnorm',Pquadfnorm);
    
    Xpointtrue = pdfnorm.transForms.normX2trueX(Xpoint);
    Xpointprior_norm = pdfnormprior.transForms.trueX2normX(Xpointtrue);
    priorprob_norm =  pdfnormprior.func(Xpointprior_norm);
    
    
%     pdfprobs_norm = zeros(size(Xx,2),1);
%     QuadFilprobs_norm = zeros(size(Xx,2),1);
%     Xpoint = zeros(1,dim);
%     for j=1:size(Xx,2)
%         Xpoint(states2keep) = [Xx(i,j),Xy(i,j)];
%         Xpoint(states2remove) = mX(states2remove);
% %         pdfprobs_norm(j) = pdfnorm.func(Xpoint);
% %         QuadFilprobs_norm(j) = mvnpdf(Xpoint,mquadfnorm',Pquadfnorm);
%         
%         pdfprobs_norm(j) = marginalize_exp_pdf_modf([Xx(i,j),Xy(i,j)],states2remove,pdfnorm,X,probs,mX,PX,GMM,'GMM_MC');
%         QuadFilprobs_norm(j) = mvnpdf([Xx(i,j),Xy(i,j)],mquadfnorm(states)',Pquadfnorm(states,states));
%     end
    pdfprobs_norm_cell{i}=pdfprobs_norm;
    QuadFilprobs_norm_cell{i}=QuadFilprobs_norm;
    priorpdfprobs_norm_cell{i}=priorprob_norm;
end

pdfprobs_norm = zeros(size(Xx));
priorpdfprobs_norm = zeros(size(Xx));
QuadFilprobs_norm = zeros(size(Xx));
for i=1:size(Xx,1)
    for j=1:size(Xx,2)
        pdfprobs_norm(i,j) = pdfprobs_norm_cell{i}(j);
        priorpdfprobs_norm(i,j) = priorpdfprobs_norm_cell{i}(j);
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


%% Get the z-probs p(z)
% [x,w] = UT_sigmapoints(mquadf,Pquadf,2);
Z=zeros(size(X,1),model.hn);
for i=1:length(probs)
   Z(i,:) = model.h(X(i,:)); 
end
[mz,Pz] = MeanCov(Z,probs/sum(probs));
Pz=Pz+model.R;

sqPz = sqrtm(Pz);
% keyboard
if model.hn==1
    Zx = linspace(mz-4*sqPz,mz+4*sqPz,100);
    Pzprobs_true = zeros(size(Zx));
    parfor i=1:1:length(Zx)
        I = integratorFuncTrueX_usingpdfnorm(pdfnormprior,@(x)mvnpdf(repmat(Zx(i),size(x,1),1),model.hvec(x),model.R),'RegTreeBoxIntegrator');
        Pzprobs_true(i) = I;
    end
    Pzk = interp1(Zx,Pzprobs_true,zk);
    figure(93)
    plot(Zx,Pzprobs_true,'b',zk,0,'k^',[zk,zk],[0,Pzk],'k--','linewidth',2,'MarkerSize',6)
    xlabel('z')
    ylabel('p(z)')
    hold off
    if saveprops.saveit==1
        saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'png')
        saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'fig')
    end
    
end
if model.hn==2
    [Zx,Zy] = meshgrid(linspace(mz(1)-4*sqPz(1,1),mz(1)+4*sqPz(1,1),19),linspace(mz(2)-4*sqPz(2,2),mz(2)+4*sqPz(2,2),19));
    Pzprobs_true = zeros(size(Zx));
    J=cell(size(Zx,1),1);
    parfor i=1:1:size(Zx,1)
        J{i}=zeros(1,size(Zx,2));
        VecFunctrue=cell(size(Zx,2),1);
        for j=1:1:size(Zx,2)
            zz=[Zx(i,j),Zy(i,j)];
            VecFunctrue{j}=@(x)mvnpdf(repmat(zz,size(x,1),1),model.hvec(x),model.R);
%             tic
%             I = integratorFuncTrueX_usingpdfnorm(pdfnormprior,@(x)mvnpdf(repmat(zz,size(x,1),1),model.hvec(x),model.R),'RegTreeBoxIntegrator');
%             toc
            %             Pzprobs_true(i,j) = I;
%             J{i}(j) = I;
        end
        tic
        J{i} = integratorVecFuncTrueX_usingpdfnorm(pdfnormprior,VecFunctrue,'RegTreeBoxIntegrator');
        toc
    end
    for i=1:1:size(Zx,1)
        for j=1:1:size(Zx,2)
            Pzprobs_true(i,j) = J{i}(j);
        end
    end
    figure(93)
    contour(Zx,Zy,Pzprobs_true,15)
    hold on
    plot(zk(1),zk(2),'k*','linewidth',2,'MarkerSize',6)
    title('Pz true contour')
    xlabel('z_1')
    ylabel('z_2')
    hold off
    if saveprops.saveit==1
        saveas(gcf,[saveprops.plotfolder,'/PZtrueContour_',nametag,'_',num2str(Tk)],'png')
        saveas(gcf,[saveprops.plotfolder,'/PZtrueContour_',nametag,'_',num2str(Tk)],'fig')
    end
    
    
    figure(94)
    surf(Zx,Zy,Pzprobs_true,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
    camlight right; lighting phong
    alpha 0.4
    hold on
    plot(zk(1),zk(2),'k*','linewidth',2,'MarkerSize',6)
    title('Pz true surf')
    xlabel('z_1')
    ylabel('z_2')
    hold off
    if saveprops.saveit==1
        saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'png')
        saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'fig')
    end

end

%%

figure(71)
contour(Xx,Xy,pdfprobs_norm,20)
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
    figure(71)
    saveas(gcf,[saveprops.plotfolder,'/NormContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(72)
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
figure(73)
contour(Xx,Xy,pdfprobs_norm,20,'r')
hold on
contour(Xx,Xy,priorpdfprobs_norm,20,'b')
grid on
box off
% h.ContourZLevel = 1;
% view([26,43])
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
    saveas(gcf,[saveprops.plotfolder,'/NormPostPriorContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormPostPriorContour_',nametag,'_',num2str(Tk)],'fig')
end
%%

figure(74)
surf(Xx,Xy,pdfprobs_norm,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
surf(Xx,Xy,priorpdfprobs_norm,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
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
    saveas(gcf,[saveprops.plotfolder,'/NormPostPriorSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormPostPriorSurf_',nametag,'_',num2str(Tk)],'fig')
end

%%
% -------------------- Quad Filter Plots------------------------------------------------------------------------
%%
%%
figure(75)
contour(Xx,Xy,QuadFilprobs_norm,20)
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
figure(76)
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