function plorprior_treepdf_4D(k,model,Xtruplot,Xtruth,priorfullnormpdf,Xmctest,Xnptk1_prior,probntk1_prior,X_pf,xfquad,Pfquad,saveprops,states)
nametag='prior';
if isempty(Xmctest)==0
    Xmctestnorm = priorfullnormpdf.transForms.trueX2normX(Xmctest);
end
Xpfnorm = priorfullnormpdf.transForms.trueX2normX(X_pf');
Xntruth = priorfullnormpdf.transForms.trueX2normX(Xtruth(k,:));
[Xq,Wq] = model.quadfunc(xfquad,Pfquad);
Xqnorm = priorfullnormpdf.transForms.trueX2normX(Xq);
[xfqnorm,Pfqnorm] = MeanCov(Xqnorm,Wq);


% figure('Name','60','Position', get(0, 'Screensize'))
% figure(60)
% if isempty(Xmctest)==0
%     plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
% end
% hold on
% plot3Dtreepatches(priorfullnormpdf.boxes)
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


figure(61)
if isempty(Xmctest)==0
    plot(Xnptk1_prior(:,1),Xnptk1_prior(:,2),'bo',Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
else
%     plot(Xnptk1_prior(:,1),Xnptk1_prior(:,2),'bo')
end
hold on
% plot_nsigellip(xfqnorm,Pfqnorm,1,'m',2)
plot2Dboxes_modf(priorfullnormpdf.boxes)
hold off
% axis([-2,2,-2,2])
grid on
xlabel('y_1')
ylabel('y_2')
% plot_paper_optim
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormTreeBox_',nametag,'_',num2str(k)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormTreeBox_',nametag,'_',num2str(k)],'fig')
end
hold off

figure(62)
if isempty(Xmctest)==0
    plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+')
end
hold on
plot_nsigellip(xfqnorm(1:2),Pfqnorm(1:2,1:2),1,'m',2)
hold off
% axis([-2,2,-2,2])
grid on
xlabel('y_1')
ylabel('y_2')
% plot_paper_optim
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormQuad_',nametag,'_',num2str(k)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormQuad_',nametag,'_',num2str(k)],'fig')
end
hold off

figure(63)
if isempty(Xmctest)==0
    plot(Xmctestnorm(:,1),Xmctestnorm(:,2),'r+',Xpfnorm(:,1),Xpfnorm(:,2),'gs')
else
%     plot(Xpfnorm(:,1),Xpfnorm(:,2),'gs')
end
% axis([-2,2,-2,2])
grid on
xlabel('y_1')
ylabel('y_2')
% plot_paper_optim
if saveprops.saveit==1
    saveas(gcf,[saveprops.plotfolder,'/NormPF_',nametag,'_',num2str(k)],'png')
    saveas(gcf,[saveprops.plotfolder,'/NormPF_',nametag,'_',num2str(k)],'fig')
end
hold off

% contour plot in norm space
[xx,yy]=meshgrid(linspace(-2,2,51));
ppnorm = zeros(size(xx));
pptruth = zeros(size(xx));
xxtruth = zeros(size(xx));
yytruth = zeros(size(xx));
ppquadtruth = zeros(size(xx));
ppquadnorm = zeros(size(xx));
NN=1e5;
intrpoints = mvurnd(-2*ones(1,2),2*ones(1,2),NN);
for i=1:size(xx,1)
    for j=1:size(xx,2)
        [i,j]
        xp=zeros(NN,model.fn);
        for r=1:NN
            xp(r,:) = [xx(i,j),yy(i,j),intrpoints(r,:)];
        end
            ssp = priorfullnormpdf.transForms.normX2trueX(xp(1,:));
            xxtruth(i,j) = ssp(1);
            yytruth(i,j) = ssp(2);
            
            s = priorfullnormpdf.normeval(xp);
            pptruth(i,j) = sum(priorfullnormpdf.transForms.normprob2trueprob(s))/NN;
            ppnorm(i,j) = sum(s)/NN;
            
            s = mvnpdf(xp,xfqnorm(:)',Pfqnorm);
            ppquadtruth(i,j) = sum(priorfullnormpdf.transForms.normprob2trueprob(s))/NN;
            ppquadnorm(i,j) = sum(s)/NN;
        
    end
end
pptruth=pptruth/16;
ppnorm=ppnorm/16;
ppquadtruth=ppquadtruth/16;
ppquadnorm=ppquadnorm/16;

% norm contour
cc = linspace(1.5*min(min(ppnorm)),max(max(ppnorm)),20);

figure(64)
contour(xx,yy,ppnorm,20)
hold on
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
figure(65)
contour(xx,yy,ppquadnorm,20)
hold on
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
figure(66)
surf(xx,yy,ppnorm,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
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
figure(67)
surf(xx,yy,ppquadnorm,'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
camlight right; lighting phong
alpha 0.4
hold on
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

%%
figure(68)
contour(xxtruth,yytruth,pptruth,20)
hold on
plot(Xtruplot(:,1),Xtruplot(:,2),'c')
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
figure(69)
contour(xxtruth,yytruth,ppquadtruth,20)
hold on
plot(Xtruplot(:,1),Xtruplot(:,2),'c')
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

