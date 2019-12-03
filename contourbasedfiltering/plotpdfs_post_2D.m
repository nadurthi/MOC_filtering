function plotpdfs_post_2D(Tk,pdfnormprior,pdfnorm,model,XtruthplotAll,X,probs,mquadf,Pquadf,mquadfprior,Pquadfprior,X_pf, w_pf,X_pfprior,w_pfprior,GMM,Xmc,Xtruth,zk,saveprops)
nametag='post';
states=[1,2];
states2keep=states;
Xn = pdfnorm.transForms.trueX2normX(X);
pn = pdfnorm.transForms.trueprob2normprob(probs);

dim = model.fn;

if isempty(Xmc)==0
    Xmcnorm = pdfnorm.transForms.trueX2normX(Xmc) ;
end
if isempty(Xtruth)==0
    Xtruthnorm = pdfnorm.transForms.trueX2normX(Xtruth) ;
end
% if isempty(mquadf)==0
%     [x,w] = UT_sigmapoints(mquadf,Pquadf,2);
%     x = pdfnorm.transForms.trueX2normX(x) ;
%     [mquadfnorm,Pquadfnorm]=MeanCov(x,w);
% end
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
Pquadf(1,2)=Pquadf(2,1);

[Xx,Xy]=meshgrid(linspace(-2,2,50),linspace(-2,2,50) );
pdfprobs_norm = zeros(size(Xx));
QuadFilprobs_norm = zeros(size(Xx));
% GMMprobs_norm = zeros(size(Xx));
for i=1:size(Xx,1)
    for j=1:size(Xx,2)
        norm_pt = [Xx(i,j),Xy(i,j)];
        pdfprobs_norm(i,j) = pdfnorm.func(norm_pt);
       
        true_pt = pdfnorm.transForms.normX2trueX(norm_pt);
        

        ptrue = mvnpdf(true_pt,mquadf(:)',Pquadf);
        QuadFilprobs_norm(i,j) = pdfnorm.transForms.trueprob2normprob(ptrue);

%         ptrue = evalGMM(GMM,true_pt);
%         GMMprobs_norm(i,j) = pdfnorm.transForms.trueprob2normprob(ptrue);
        
    end
end
normmaxprob = max([max(max(pdfprobs_norm)), max(max(QuadFilprobs_norm))]);
normminprob = min([min(min(pdfprobs_norm)), min(min(QuadFilprobs_norm))]);
normproblinspace = linspace(1.1*normminprob,0.9*normmaxprob,100);

% get the true points and their probs
pdfprobs_true = zeros(size(Xx));
priorpdfprobs_true = zeros(size(Xx));
QuadFilprobs_true = zeros(size(Xx));
% GMMprobs_true = zeros(size(Xx));
Xx_true = zeros(size(Xx));
Xy_true = zeros(size(Xx));

[a,b]=size(Xx);

kspfpdfnormpost = mvksdensity(X_pfnorm(states2keep,:)',[reshape(Xx,a*b,1),reshape(Xy,a*b,1)],...
	'Bandwidth',bwpost(states2keep),...
	'Kernel','normpdf');

kspfpdfnormprior = mvksdensity(X_pfpriornorm(states2keep,:)',[reshape(Xx,a*b,1),reshape(Xy,a*b,1)],...
	'Bandwidth',bwprior(states2keep),...
	'Kernel','normpdf');

kspfpdfnormpost=reshape(kspfpdfnormpost,a,b);
kspfpdfnormprior=reshape(kspfpdfnormprior,a,b);

for i=1:size(Xx,1)
    for j=1:size(Xx,2)
        pttrue=pdfnorm.transForms.normX2trueX([Xx(i,j),Xy(i,j)]);
        Xx_true(i,j) = pttrue(1);
        Xy_true(i,j) = pttrue(2);
        pdfprobs_true(i,j) = pdfnorm.transForms.normprob2trueprob( pdfprobs_norm(i,j) );
        
        ptnorm_prior = pdfnormprior.transForms.trueX2normX(pttrue);
        prior_prob_norm = pdfnormprior.func(ptnorm_prior);
        priorpdfprobs_true(i,j) = pdfnormprior.transForms.normprob2trueprob( prior_prob_norm );
        
        QuadFilprobs_true(i,j) = mvnpdf([Xx_true(i,j),Xy_true(i,j)],mquadf(:)',Pquadf);
%         GMMprobs_true(i,j) = evalGMM(GMM,[Xx_true(i,j),Xy_true(i,j)]);
        
    end
end



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
    Zx = linspace(mz-3*sqPz,mz+3*sqPz,100);
%     Pzprobs_true = zeros(size(Zx));
%     for i=1:1:length(Zx)
%         I = integratorFuncTrueX_usingpdfnorm(pdfnormprior,@(x)mvnpdf(repmat(Zx(i),size(x,1),1),model.h(x),model.R),'RegTreeBoxIntegrator');
%         Pzprobs_true(i) = I;
%     end
    Pzprobs_true = zeros(size(Zx));
%     J=cell(size(Zx,1),1);
    LB=-1.1*ones(1,dim);
    UB=1.1*ones(1,dim);
    voll = prod(UB-LB);
    [XX,WW]=GLgn_pts(LB,UB,9);
    %         keyboard
    for i=1:1:length(Zx)
        i

            zz=Zx(i);
            pnorm = pdfnormprior.func(XX);
            Xtr=pdfnormprior.transForms.normX2trueX(XX);
            pz=mvnpdf(repmat(zz,size(Xtr,1),1),model.hvec(Xtr),model.R);
            Pzprobs_true(i)=sum(WW.*pnorm.*pz);


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
    [Zx,Zy] = meshgrid(linspace(mz(1)-3*sqPz(1,1),mz(1)+3*sqPz(1,1),25),linspace(mz(2)-3*sqPz(2,2),mz(2)+3*sqPz(2,2),25));
    Pzprobs_true = zeros(size(Zx));
    J=cell(size(Zx,1),1);
    parfor i=1:1:size(Zx,1)
        J{i}=zeros(1,size(Zx,2));
        for j=1:1:size(Zx,2)
            zz=[Zx(i,j),Zy(i,j)];
            I = integratorFuncTrueX_usingpdfnorm(pdfnormprior,@(x)mvnpdf(repmat(zz,size(x,1),1),model.h(x),model.R),'RegTreeBoxIntegrator');
%             Pzprobs_true(i,j) = I;
            J{i}(j) = I;
        end
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
% keyboard
figure(41)
contour(Xx,Xy,pdfprobs_norm,50)
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
title('Post Reg Contour norm')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(42)
surf(Xx,Xy,pdfprobs_norm,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
plot3(Xn(:,1),Xn(:,2),pn,'bo')

if isempty(Xmc)==0
    plot(Xmcnorm(:,1),Xmcnorm(:,2),'r.')
end
if isempty(Xtruthnorm)==0
    plot(Xtruthnorm(:,1),Xtruthnorm(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Post Reg Surf norm')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormSurf_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(43)
[~,h]=contour(Xx_true,Xy_true,pdfprobs_true ,15);
grid on
box off
% h.ContourZLevel = 1;
% view([26,43])
hold on
if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Post Reg Contour true')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/TrueContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/TrueContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(44)
surf(Xx_true,Xy_true,pdfprobs_true,'FaceColor','g','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.7
hold on
plot3(X(:,1),X(:,2),probs,'bo')

if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Post Reg Surf true')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/TrueSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/TrueSurf_',nametag,'_',num2str(Tk)],'fig')
end

%%
% -------------------- Quad Filter Plots------------------------------------------------------------------------
%%
%%
figure(45)
contour(Xx,Xy,QuadFilprobs_norm,50)
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
title('Post Quad Contour norm')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(46)
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
title('Post Quad Surf norm')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilNormSurf_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(47)
[~,h]=contour(Xx_true,Xy_true,QuadFilprobs_true ,15);
grid on
box off
% h.ContourZLevel = 1;
% view([26,43])
hold on
if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Post Quad Contour true')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueContour_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueContour_',nametag,'_',num2str(Tk)],'fig')
end
%%
figure(48)
surf(Xx_true,Xy_true,QuadFilprobs_true,'FaceColor','r','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.7
hold on
% plot3(X(:,1),X(:,2),probs,'bo')

if isempty(Xmc)==0
    plot(Xmc(:,1),Xmc(:,2),'r.')
end
if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

% title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
xlabel('x_1')
ylabel('x_2')
title('Post Quad Surf true')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/QuadFilTrueSurf_',nametag,'_',num2str(Tk)],'fig')
end

%% PLot pdfs on full true trajectory
figure(49)
plot(XtruthplotAll(:,1),XtruthplotAll(:,2),'k--')
hold on
box off
grid on

[~,h]=contour(Xx_true,Xy_true,priorpdfprobs_true ,15,'LineColor','b');
[~,h]=contour(Xx_true,Xy_true,pdfprobs_true ,15,'LineColor','r');

xlabel('x_1')
ylabel('x_2')
title('Post and Prior Contour')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/FullTrajTrueCont_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/FullTrajTrueCont_',nametag,'_',num2str(Tk)],'fig')
end

figure(50)
plot(XtruthplotAll(:,1),XtruthplotAll(:,2),'k--')
hold on
box off
grid on

surf(Xx_true,Xy_true,priorpdfprobs_true,'FaceColor','b','EdgeColor','none','FaceAlpha',0.5);
camlight right; lighting phong
alpha 0.5

surf(Xx_true,Xy_true,pdfprobs_true,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
camlight right; lighting phong
alpha 0.5

xlabel('x_1')
ylabel('x_2')
title('Post and Prior Surf')
axis equal
axis square
hold off
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/FullTrajTrueSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/FullTrajTrueSurf_',nametag,'_',num2str(Tk)],'fig')
end

%% PLot pdfs both prior and post
figure(51)
% plot(XtruthplotAll(:,1),XtruthplotAll(:,2),'k--')


[~,h]=contour(Xx_true,Xy_true,priorpdfprobs_true ,15,'LineColor','b');
hold on
box off
grid on

[~,h]=contour(Xx_true,Xy_true,pdfprobs_true ,15,'LineColor','r');

if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

xlabel('x_1')
ylabel('x_2')
title('Post and Prior Contour')
axis equal
axis square
hold off


if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/BothTrueCont_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/BothTrueCont_',nametag,'_',num2str(Tk)],'fig')
end

figure(52)
% plot(XtruthplotAll(:,1),XtruthplotAll(:,2),'k--')


surf(Xx_true,Xy_true,priorpdfprobs_true,'FaceColor','b','EdgeColor','none','FaceAlpha',0.5);
camlight right; lighting phong
alpha 0.5

hold on
box off
grid on

surf(Xx_true,Xy_true,pdfprobs_true,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
camlight right; lighting phong
alpha 0.5

if isempty(Xtruth)==0
    plot(Xtruth(:,1),Xtruth(:,2),'k*','linewidth',2)
end

xlabel('x_1')
ylabel('x_2')
title('Post and Prior Surf')
axis equal
axis square
hold off

if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/BothTrueSurf_',nametag,'_',num2str(Tk)],'png')
    saveas(gcf,[saveprops.plotfolder,'/BothTrueSurf_',nametag,'_',num2str(Tk)],'fig')
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

