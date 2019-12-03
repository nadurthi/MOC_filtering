function plotpdfs_post_6D(states,Tk,pdfnormprior,pdfnorm,model,X,probs,mquadf,Pquadf,mquadfprior,Pquadfprior,X_pf, w_pf,X_pfprior,w_pfprior,Xmc,Xtruth,zk,saveprops,doPZ)
tagcoordstr=10*states(1)+states(2);
tagcoordstr = num2str(tagcoordstr);
nametag=strcat([tagcoordstr,'_prior']);

%%
dim=size(X,2);
states2keep = states;
states2remove = 1:dim;
states2remove(states2keep)=[];


% Xn = pdfnorm.transForms.trueX2normX(X);
% pn = pdfnorm.transForms.trueprob2normprob(probs);
% GMMfitter = GMMFitDataSet(Xn,pn);
% GMM = GMMfitter.FitGMM_kmeans_optimwt(3);
%% resample pf to make equal weightd
Npf=length(w_pf);
I = myresampling(w_pf);
I = round(I);
X_old=X_pf;
for j = 1 : Npf
    X_pf(:,j) = X_old(:,I(j));
end
w_pf = ones(1,Npf)/Npf;

I = myresampling(w_pfprior);
I = round(I);
X_old=X_pfprior;
for j = 1 : Npf
    X_pfprior(:,j) = X_old(:,I(j));
end
w_pfprior = ones(1,Npf)/Npf;

X_pfnorm=zeros(size(X_pf));
X_pfpriornorm=zeros(size(X_pf));
for i=1:Npf
 X_pfnorm(:,i)=pdfnorm.transForms.trueX2normX(X_pf(:,i)');
 X_pfpriornorm(:,i)=pdfnorm.transForms.trueX2normX(X_pfprior(:,i)');
end
bwpost=zeros(1,dim);
bwprior=zeros(1,dim);

Npf = length(w_pf);
for i=1:dim
    s=std(X_pfnorm(:,i));
    bwpost(i) = s*(4/((dim+2)*Npf))^(1/(dim+4));
    
    s=std(X_pfpriornorm(:,i));
    bwprior(i) = s*(4/((dim+2)*Npf))^(1/(dim+4));
    
end
%%

if isempty(Xmc)==0
    Xmcnorm = pdfnorm.transForms.trueX2normX(Xmc) ;
end
if isempty(Xtruth)==0
    Xtruthnorm = pdfnorm.transForms.trueX2normX(Xtruth) ;
end
if isempty(mquadf)==0
    [x,w] = UT_sigmapoints(mquadf,Pquadf,2);
    x = pdfnorm.transForms.trueX2normX(x) ;
    [mquadfnormpost,Pquadfnormpost]=MeanCov(x,w);
    mquadfnormpost=mquadfnormpost(:);
    
    [x,w] = UT_sigmapoints(mquadfprior,Pquadfprior,2);
    x = pdfnorm.transForms.trueX2normX(x) ;
    [mquadfnormprior,Pquadfnormprior]=MeanCov(x,w);
    mquadfnormprior=mquadfnormprior(:);
    
end

% [mXn,PXn] = MeanCov(Xn,pn/sum(pn));

[Xx,Xy]=meshgrid(linspace(-2,2,25),linspace(-2,2,25) );
% [Xx,Xy]=meshgrid(linspace(-1.4,1.4,35),linspace(-1.4,1.4,35) );
%%
QuadFilprobspost_norm = zeros(size(Xx));
QuadFilprobsprior_norm = zeros(size(Xx));

for i=1:size(Xx,1)
    for j=1:size(Xx,2)
        QuadFilprobspost_norm(i,j)=mvnpdf([Xx(i,j),Xy(i,j)],mquadfnormpost(states2keep)',Pquadfnormpost(states2keep,states2keep));
        QuadFilprobsprior_norm(i,j)=mvnpdf([Xx(i,j),Xy(i,j)],mquadfnormprior(states2keep)',Pquadfnormprior(states2keep,states2keep));
    end  
end

remdim=length(states2remove);
remdomainLB=-1.1*ones(remdim,1);
remdomainUB=1.1*ones(remdim,1);
volremdom = prod(remdomainUB-remdomainLB);
[Xstates2remove,wrem]=GLgn_pts(remdomainLB,remdomainUB,9);


pdfprobsprior_norm_cell = cell(size(Xx,1),1);
pdfprobspost_norm_cell = cell(size(Xx,1),1);

% priorpdfprobs_norm_cell = cell(size(Xx,1),1);
% QuadFilprobs_norm_cell = cell(size(Xx,1),1);

[a,b]=size(Xx);

kspfpdfnormpost = mvksdensity(X_pfnorm(states2keep,:)',[reshape(Xx,a*b,1),reshape(Xy,a*b,1)],...
	'Bandwidth',bwpost(states2keep),...
	'Kernel','normpdf');

kspfpdfnormprior = mvksdensity(X_pfpriornorm(states2keep,:)',[reshape(Xx,a*b,1),reshape(Xy,a*b,1)],...
	'Bandwidth',bwprior(states2keep),...
	'Kernel','normpdf');

kspfpdfnormpost=reshape(kspfpdfnormpost,a,b);
kspfpdfnormprior=reshape(kspfpdfnormprior,a,b);

parfor i=1:size(Xx,1)
    i
    pmargnormpost = zeros(size(Xx,2),1);
    pmargnormprior = zeros(size(Xx,2),1);
     
    for j=1:size(Xx,2)
        Xpoint = zeros(size(Xstates2remove,1),dim);
        Xpoint(:,states2keep) = repmat([Xx(i,j),Xy(i,j)],size(Xstates2remove,1),1);
        Xpoint(:,states2remove) = Xstates2remove;
        
        ff= pdfnorm.func(Xpoint);
        pmargnormpost(j) = wrem(:)'*ff(:);  
        
        xtrpost = pdfnorm.transForms.normX2trueX(Xpoint) ;
        xnormprior = pdfnormprior.transForms.trueX2normX(xtrpost) ;
        
        ffpriornorm = pdfnormprior.func(xnormprior);
        ffpriortrue = pdfnormprior.transForms.normprob2trueprob(ffpriornorm) ;
        ffpostnorm = pdfnorm.transForms.trueprob2normprob(ffpriortrue) ;
        
        pmargnormprior(j) = wrem(:)'*ffpostnorm(:); 
        
%         if pmargnormpost(j)>10
%             keyboard
%         end
        
    end
    pmargnormpost=volremdom*pmargnormpost;
    pmargnormprior=volremdom*pmargnormprior;

    pdfprobspost_norm_cell{i}=pmargnormpost;
    pdfprobsprior_norm_cell{i}=pmargnormprior;
end
% keyboard

pdfprobspost_norm = zeros(size(Xx));
pdfprobsprior_norm = zeros(size(Xx));

for i=1:size(Xx,1)
    pdfprobspost_norm(i,:) = pdfprobspost_norm_cell{i};
    pdfprobsprior_norm(i,:) = pdfprobsprior_norm_cell{i};
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

if doPZ==1
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
        Zx = linspace(mz-3*sqPz,mz+3*sqPz,100);
        Pzprobs_true = zeros(size(Zx));
        for i=1:1:length(Zx)
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
        disp('-----------------doing PZ-----------------------')
        [Zx,Zy] = meshgrid(linspace(mz(1)-4*sqPz(1,1),mz(1)+4*sqPz(1,1),19),linspace(mz(2)-4*sqPz(2,2),mz(2)+4*sqPz(2,2),19));
        Pzprobs_true = zeros(size(Zx));
        J=cell(size(Zx,1),1);
        LB=-1.1*ones(1,dim);
        UB=1.1*ones(1,dim);
        voll = prod(UB-LB);
        [XX,WW]=GLgn_pts(LB,UB,4);
%         keyboard
        parfor i=1:1:size(Zx,1)
            i
            J{i}=zeros(1,size(Zx,2));
            for j=1:1:size(Zx,2)
                zz=[Zx(i,j),Zy(i,j)];
                pnorm = pdfnormprior.func(XX);
                Xtr=pdfnormprior.transForms.normX2trueX(XX);
                pz=mvnpdf(repmat(zz,size(Xtr,1),1),model.hvec(Xtr),model.R);
                J{i}(j)=sum(WW.*pnorm.*pz);
            end

        end
        for i=1:1:size(Zx,1)
            Pzprobs_true(i,:) = J{i};
        end
        figure(93)
        contour(Zx,Zy,Pzprobs_true,15,'linewidth',1)
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
        xlabel(model.hstates{1})
        ylabel(model.hstates{2})
        hold off
        if saveprops.saveit==1
            saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'png')
            saveas(gcf,[saveprops.plotfolder,'/PZtrueSurf_',nametag,'_',num2str(Tk)],'fig')
        end
        
    end
    
end
%%
mx=max(max(pdfprobspost_norm));
mn=min(min(pdfprobspost_norm));
ctrs = linspace(mn,mx,20);

figure(71)
contour(Xx,Xy,pdfprobspost_norm,ctrs,'linewidth',1)
grid on
box off
hold on
if isempty(Xmc)==0
    plot(Xmcnorm(:,states(1)),Xmcnorm(:,states(2)),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,states(1)),Xtruthnorm(:,states(2)),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel(model.fstates{states(1)})
ylabel(model.fstates{states(2)})
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(72)
surf(Xx,Xy,pdfprobspost_norm,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
% plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,states(1)),Xmcnorm(:,states(2)),'r.')
end
if isempty(Xtruthnorm)==0
    plot(Xtruthnorm(:,states(1)),Xtruthnorm(:,states(2)),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel(model.fstates{states(1)})
ylabel(model.fstates{states(2)})
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormSurf_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(73)
contour(Xx,Xy,pdfprobspost_norm,ctrs,'r','linewidth',1)
hold on
contour(Xx,Xy,pdfprobsprior_norm,ctrs,'b','linewidth',1)
grid on
box off
% h.ContourZLevel = 1;
% view([26,43])
hold on
if isempty(Xmc)==0
    plot(Xmcnorm(:,states(1)),Xmcnorm(:,states(2)),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,states(1)),Xtruthnorm(:,states(2)),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel(model.fstates{states(1)})
ylabel(model.fstates{states(2)})
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormPostPriorContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormPostPriorContour_',nametag,'_',num2str(Tk)],'fig')
end
%%

figure(74)
surf(Xx,Xy,pdfprobspost_norm,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
surf(Xx,Xy,pdfprobsprior_norm,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on

if isempty(Xmc)==0
    plot(Xmcnorm(:,states(1)),Xmcnorm(:,states(2)),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,states(1)),Xtruthnorm(:,states(2)),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel(model.fstates{states(1)})
ylabel(model.fstates{states(2)})
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
contour(Xx,Xy,QuadFilprobspost_norm,ctrs,'r','linewidth',1)
hold on
contour(Xx,Xy,QuadFilprobsprior_norm,ctrs,'b','linewidth',1)

grid on
box off
hold on
if isempty(Xmc)==0
    plot(Xmcnorm(:,states(1)),Xmcnorm(:,states(2)),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,states(1)),Xtruthnorm(:,states(2)),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel(model.fstates{states(1)})
ylabel(model.fstates{states(2)})
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormPostPriorContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormPostPriorContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(76)
surf(Xx,Xy,QuadFilprobspost_norm,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on

surf(Xx,Xy,QuadFilprobsprior_norm,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4

% plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,states(1)),Xmcnorm(:,states(2)),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,states(1)),Xtruthnorm(:,states(2)),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel(model.fstates{states(1)})
ylabel(model.fstates{states(2)})
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormPostPriorSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormPostPriorSurf_',nametag,'_',num2str(Tk)],'fig')
end
%% -------------------- KS-density pf ------------------------------------------------------------------------



figure(77)
contour(Xx,Xy,kspfpdfnormpost,20,'r','linewidth',1)
hold on
contour(Xx,Xy,kspfpdfnormprior,20,'b','linewidth',1)

grid on
box off
hold on
if isempty(Xmc)==0
    plot(Xmcnorm(:,states(1)),Xmcnorm(:,states(2)),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,states(1)),Xtruthnorm(:,states(2)),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel(model.fstates{states(1)})
ylabel(model.fstates{states(2)})
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/kspfFilNormPostPriorContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/kspfFilNormPostPriorContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(78)
surf(Xx,Xy,kspfpdfnormpost,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on

surf(Xx,Xy,kspfpdfnormprior,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4

% plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,states(1)),Xmcnorm(:,states(2)),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruthnorm(:,states(1)),Xtruthnorm(:,states(2)),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel(model.fstates{states(1)})
ylabel(model.fstates{states(2)})
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/kspfFilNormPostPriorSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/kspfFilNormPostPriorSurf_',nametag,'_',num2str(Tk)],'fig')
end

% keyboard
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