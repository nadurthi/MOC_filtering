function plotpost_treepdf_6D(k,model,Xtruplot,Xtruth,priorfullnormpdf,postfullnormpdf,Xmctest,Xnptk1_prior,probntk1_prior,Xnptk1_post,probntk1_post,X_pf,w_pf,xfquadprior,Pfquadprior,xfquad,Pfquad,saveprops,plotstates)
nametag='post';
if isempty(Xmctest)==0
    Xmctestnorm = priorfullnormpdf.transForms.trueX2normX(Xmctest);
end

% X_pf,w_pf
X_pfnew = zeros(size(X_pf));
I = myresampling(w_pf);
I = round(I);
for j = 1 : length(w_pf)
    X_pfnew(:,j) = X_pf(:,I(j));
end
% reset weights
w_new = ones(1,length(w_pf))/length(w_pf);
X_pf = X_pfnew;

Xpfnorm = priorfullnormpdf.transForms.trueX2normX(X_pf');
Xntruth = priorfullnormpdf.transForms.trueX2normX(Xtruth(k,:));
[Xq,Wq] = model.quadfunc(xfquad,Pfquad);
Xqnorm = priorfullnormpdf.transForms.trueX2normX(Xq);
[xfqnorm,Pfqnorm] = MeanCov(Xqnorm,Wq);

[Xqprior,Wqprior] = model.quadfunc(xfquadprior,Pfquadprior);
Xqpriornorm = priorfullnormpdf.transForms.trueX2normX(Xqprior);
[xfqnormprior,Pfqnormprior] = MeanCov(Xqpriornorm,Wqprior);

% figure('Name','60','Position', get(0, 'Screensize'))
% figure(50)
% if isempty(Xmctest)==0
%     plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
% end
% hold on
% plot3Dtreepatches_col(priorfullnormpdf.boxes,'b')
% plot3Dtreepatches_col(postfullnormpdf.boxes,'r')
% plot(Xpfnorm(:,1),Xpfnorm(:,2),'gs')
% 
% hold off
% grid on;
% % axis([-2,2,-2,2])
% xlabel('y_1')
% ylabel('y_2')
% zlabel('p(y)')
% % plot_paper_optim
% if saveprops.saveit==1
%     saveas(gcf,[saveprops.plotfolder,'/NormTreeHeight_',nametag,'_',num2str(k)],'png')
%     saveas(gcf,[saveprops.plotfolder,'/NormTreeHeight_',nametag,'_',num2str(k)],'fig')
% end



% figure(51)
% if isempty(Xmctest)==0
%     plot(Xnptk1_prior(:,1),Xnptk1_prior(:,2),'bo',Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
% else
%     plot(Xnptk1_prior(:,1),Xnptk1_prior(:,2),'bo')
% end
% hold on
% % plot_nsigellip(xfqnorm,Pfqnorm,1,'m',2)
% plot2Dboxes_modf(priorfullnormpdf.boxes)
% hold off
% % axis([-2,2,-2,2])
% grid on
% xlabel('y_1')
% ylabel('y_2')
% % plot_paper_optim
% if saveprops.saveit==1
%     saveas(gcf,[saveprops.plotfolder,'/NormTreeBox_',nametag,'_',num2str(k)],'png')
%     saveas(gcf,[saveprops.plotfolder,'/NormTreeBox_',nametag,'_',num2str(k)],'fig')
% end
% hold off


% 
% figure(52)
% if isempty(Xmctest)==0
%     plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
% end
% hold on
% plot_nsigellip(xfqnormprior,Pfqnormprior,1,'b',2)
% plot_nsigellip(xfqnorm,Pfqnorm,1,'r',2)
% hold off
% % axis([-2,2,-2,2])
% grid on
% xlabel('y_1')
% ylabel('y_2')
% % plot_paper_optim
% if saveprops.saveit==1
%     saveas(gcf,[saveprops.plotfolder,'/NormQuad_',nametag,'_',num2str(k)],'png')
%     saveas(gcf,[saveprops.plotfolder,'/NormQuad_',nametag,'_',num2str(k)],'fig')
% end
% hold off

% figure(53)
% if isempty(Xmctest)==0
%     plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+',Xpfnorm(:,1),Xpfnorm(:,2),'gs')
% else
%     plot(Xpfnorm(:,1),Xpfnorm(:,2),'gs')
% end
% % axis([-2,2,-2,2])
% grid on
% xlabel('y_1')
% ylabel('y_2')
% % plot_paper_optim
% if saveprops.saveit==1
%     saveas(gcf,[saveprops.plotfolder,'/NormPF_',nametag,'_',num2str(k)],'png')
%     saveas(gcf,[saveprops.plotfolder,'/NormPF_',nametag,'_',num2str(k)],'fig')
% end
% hold off

% contour plot in norm space
[xx,yy]=meshgrid([linspace(-2,-0.7,10),linspace(-0.68,0.68,101),linspace(0.7,2,10)]);
ppnorm = zeros(size(xx));
pptruth = zeros(size(xx));
ppnormpost = zeros(size(xx));
pptruthpost = zeros(size(xx));

xxtruth = zeros(size(xx));
yytruth = zeros(size(xx));

ppquadtruth = zeros(size(xx));
ppquadnorm = zeros(size(xx));
ppquadtruthpost = zeros(size(xx));
ppquadnormpost = zeros(size(xx));
for i=1:size(xx,1)
    for j=1:size(xx,2)
        
        ss = priorfullnormpdf.transForms.normX2trueX([xx(i,j),yy(i,j)]);
        xxtruth(i,j) = ss(1);
        yytruth(i,j) = ss(2);
        ppnorm(i,j) = priorfullnormpdf.normeval([xx(i,j),yy(i,j)]);
        pptruth(i,j) = priorfullnormpdf.transForms.normprob2trueprob(ppnorm(i,j));
        
        ppnormpost(i,j) = postfullnormpdf.normeval([xx(i,j),yy(i,j)]);
        pptruthpost(i,j) = postfullnormpdf.transForms.normprob2trueprob(ppnormpost(i,j));
        
        ppquadnorm(i,j) = mvnpdf([xx(i,j),yy(i,j)],xfqnormprior(:)',Pfqnormprior);
        ppquadtruth(i,j) = priorfullnormpdf.transForms.normprob2trueprob(ppquadnorm(i,j));
        
        ppquadnormpost(i,j) = mvnpdf([xx(i,j),yy(i,j)],xfqnorm(:)',Pfqnorm);
        ppquadtruthpost(i,j) = priorfullnormpdf.transForms.normprob2trueprob(ppquadnormpost(i,j));
    end
end

% norm contour
cc = linspace(1.5*min(min(ppnorm)),max(max(ppnorm)),20);

figure(54)
contour(xx,yy,ppnorm,20,'b')
hold on
contour(xx,yy,ppnormpost,20,'r')
% plot(Xpfnorm(:,1),Xpfnorm(:,2),'gs')
plot(Xntruth(:,1),Xntruth(:,2),'kx','linewidth',2,'MarkerSize',9)

if isempty(Xmctest)==0
    plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
end
% axis([-2,2,-2,2])
xlabel('y_1')
ylabel('y_2')
% plot_paper_optim

if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormTreeContour_',nametag,'_',num2str(k)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormTreeContour_',nametag,'_',num2str(k)],'fig')
end
hold off

% quad norm
figure(55)
contour(xx,yy,ppquadnorm,20,'b')
hold on
contour(xx,yy,ppquadnormpost,20,'r')
% plot(Xpfnorm(:,1),Xpfnorm(:,2),'gs')
plot(Xntruth(:,1),Xntruth(:,2),'kx','linewidth',2,'MarkerSize',9)
if isempty(Xmctest)==0
    plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
end
% axis([-2,2,-2,2])
xlabel('y_1')
ylabel('y_2')
% plot_paper_optim
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormQuadContour_',nametag,'_',num2str(k)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormQuadContour_',nametag,'_',num2str(k)],'fig')
end
hold off

%% surf plots
figure(56)
surf(xx,yy,ppnorm,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
surf(xx,yy,ppnormpost,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
% plot(Xpfnorm(:,1),Xpfnorm(:,2),'gs')
histogram2(Xpfnorm(:,1),Xpfnorm(:,2),'Normalization','pdf','FaceColor','g','FaceAlpha',0.4)
plot(Xntruth(:,1),Xntruth(:,2),'kx','linewidth',2,'MarkerSize',9)
if isempty(Xmctest)==0
    plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
end
% axis([-2,2,-2,2])
xlabel('y_1')
ylabel('y_2')
ylabel('p(y)')
% plot_paper_optim
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormTreeSurf_',nametag,'_',num2str(k)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormTreeSurf_',nametag,'_',num2str(k)],'fig')
end
hold off

% quad norm
figure(57)
surf(xx,yy,ppquadnorm,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
surf(xx,yy,ppquadnormpost,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
% plot(Xpfnorm(:,1),Xpfnorm(:,2),'gs')
histogram2(Xpfnorm(:,1),Xpfnorm(:,2),'Normalization','pdf','FaceColor','g','FaceAlpha',0.4)
plot(Xntruth(:,1),Xntruth(:,2),'kx','linewidth',2,'MarkerSize',9)
if isempty(Xmctest)==0
    plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
end
% axis([-2,2,-2,2])
xlabel('y_1')
ylabel('y_2')
ylabel('p(y)')
% plot_paper_optim
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormQuadSurf_',nametag,'_',num2str(k)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormQuadSurf_',nametag,'_',num2str(k)],'fig')
end
hold off

%% true contours
Xpf=X_pf';
figure(58)
contour(xxtruth,yytruth,pptruth,20,'b')
hold on
contour(xxtruth,yytruth,pptruthpost,20,'r')
plot(Xtruplot(:,1),Xtruplot(:,2),'c')
plot(Xpf(:,1),Xpf(:,2),'gs')
if isempty(Xmctest)==0
    plot(Xmctest(:,1),Xmctest(:,2),'r+')
end
% axis([-8,8,-25,25])
xlabel('x_1')
ylabel('x_2')
ylabel('p(x)')
% plot_paper_optim
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/TruthTreeContour_',nametag,'_',num2str(k)],'png')
    saveas(gcf,[saveprops.plotfolder,'/TruthTreeContour_',nametag,'_',num2str(k)],'fig')
end
hold off

% quad norm
figure(59)
contour(xxtruth,yytruth,ppquadtruth,20,'b')
hold on
contour(xxtruth,yytruth,ppquadtruthpost,20,'r')
plot(Xtruplot(:,1),Xtruplot(:,2),'c')
plot(Xpf(:,1),Xpf(:,2),'gs')
if isempty(Xmctest)==0
    plot(Xmctest(:,1),Xmctest(:,2),'r+')
end
% axis([-8,8,-25,25])
xlabel('x_1')
ylabel('x_2')
ylabel('p(x)')
% plot_paper_optim
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/TruthQuadContour_',nametag,'_',num2str(k)],'png')
    saveas(gcf,[saveprops.plotfolder,'/TruthQuadContour_',nametag,'_',num2str(k)],'fig')
end
hold off