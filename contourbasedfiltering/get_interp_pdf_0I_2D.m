function pdfnorm = get_interp_pdf_0I_2D(X,probs,mquad,Pquad,Nm,Tk,Xmctest,Xtruth,plotsconf)
%%

dim =size(X,2);
%
% [m,P]=MeanCov(X,probs/sum(probs));
c=1;
m=mquad(:);
P=Pquad/c;
% P=eye(dim);
%
N=size(X,1);

Z=zeros(size(X));
if isempty(Xmctest)==0
    Zmctest = zeros(size(Xmctest));
    Nmctest = size(Xmctest,1);
end
if isempty(Xtruth)==0
    Ztruth = zeros(size(Xtruth));
    Ntruth = size(Xtruth,1);
end

Psqrt=sqrtm(P);
Psqrt_inv=inv(sqrtm(P));
for i=1:N
    Z(i,:)=Psqrt_inv*(X(i,:)-m')';
end
if isempty(Xmctest)==0
    for i=1:Nmctest
        Zmctest(i,:)=Psqrt_inv*(Xmctest(i,:)-m')';
    end
end
if isempty(Xtruth)==0
    for i=1:Ntruth
        Ztruth(i,:)=Psqrt_inv*(Xtruth(i,:)-m')';
    end
end

pn=probs*det(Psqrt);

% remove points outside 10-sigma
ind = sqrt(sum(Z.^2,2))<10*sqrt(c);
Z=Z(ind,:);
pn=pn(ind);
N=size(Z,1);

% do hyopercube scalling
mn = min(Z,[],1);
Xn=zeros(size(Z));

for i=1:N
    Xn(i,:)=Z(i,:)-mn;
end
if isempty(Xmctest)==0
    Xnmctest=zeros(size(Zmctest));
    for i=1:Nmctest
        Xnmctest(i,:)=Zmctest(i,:)-mn;
    end
end
if isempty(Xtruth)==0
    Xntruth=zeros(size(Ztruth));
    for i=1:Ntruth
        Xntruth(i,:)=Ztruth(i,:)-mn;
    end
end

mx = max(Xn,[],1);
Atransf=diag(2./mx);
mulin = -2*mn(:)./mx(:)-1;
for i=1:N
    Xn(i,:)=Atransf*Xn(i,:)'-1;
end
if isempty(Xmctest)==0
    for i=1:Nmctest
        Xnmctest(i,:)=Atransf*Xnmctest(i,:)'-1;
    end
end
if isempty(Xtruth)==0
    for i=1:Ntruth
        Xntruth(i,:)=Atransf*Xntruth(i,:)'-1;
    end
end


detAtransf = det(Atransf);

% get full prob transformation
pn=pn/det(Atransf);
ind=pn<1e-70;
pn(ind)=1e-70;

% sort the points from highest to lowets prob
[~,ind]=sort(pn(:),1,'descend');
pn=pn(ind);
Xn=Xn(ind,:);

logpn = log(pn);


%% plottinmg

figure(33)
plot3(Xn(:,1),Xn(:,2),log(pn),'ro')
title(['time step = ',num2str(Tk)])
% plot3(Y(:,1),Y(:,2),log(pn),'ro')
%% fitting poly to log of probas
% keyboard

Pf=Basis_polyND(dim,Nm);
% now interpolate the polynomials to Xn and logpn
% Nleastfit = 3*length(Pf);

A=zeros(N,length(Pf));
tic
for r=1:1:N
    A(r,:) = evaluate_MatrixOfPolys(Pf,Xn(r,:));
end
lam = A\log(pn);
disp('norm fit')
[rank(A),max(abs(A*lam-log(pn))),min(lam),max(lam)]
toc
% evaluate_polyND(Pf{65},Xn(11,:))
% evaluate_polyND_3(Pf{65},Xn(11,:))




%% reguralization points
rad=max(sqrt(sum(Xn.^2,2)));
Xbnd=[]
% Xbnd=[Xbnd;2.5*gen_uniform_grid(8,dim)];
% Xbnd=GH_points(zeros(dim,1),eye(dim),6);
% Xbnd=3*Xbnd/max(max(Xbnd));
% Xbnd = [Xbnd;3*(rand(3500,dim)*2-1)];
[Xbnd1,~] = GLgn_pts(-2.5*ones(1,dim),2.5*ones(1,dim),10);
% [Xbnd2,~] = GLgn_pts(-2*ones(1,dim),2*ones(1,dim),8);
% Xbnd2=2.5*(rand(200,dim)*2-1);
Xbnd2=2.5*gen_uniform_grid(13,dim);
Xbnd = [Xbnd1;Xbnd2];
% Xbnd=1*Xbnd(sqrt(sum(Xbnd.^2,2))>1.5*rad,:);
[size(Xbnd),length(lam)]
removeind=[];
for i=1:size(Xbnd,1)
    Idx = knnsearch(Xn,Xbnd(i,:),'K',5);
    x = Xn(Idx,:);
    mr = mean(x,1);
    r  = max(sqrt(sum((x - repmat(mr,size(x,1),1)).^2)));
%     if norm(Xbnd(i,:)-mr)<0.5
    if norm(Xbnd(i,:)-mr)<2*r
        removeind=[removeind,i];
    end
%     end
end
size(Xbnd,1)
Xbnd(removeind,:)=[];
figure(35)
plot3(Xn(:,1),Xn(:,2),log(pn),'ro')
hold on
plot3(Xbnd(:,1),Xbnd(:,2),-ones(size(Xbnd,1),1),'b+')
hold off

%
M=size(Xbnd,1);
Dineq = zeros(M,length(Pf));
for r=1:1:M
    Dineq(r,:) = evaluate_MatrixOfPolys(Pf,Xbnd(r,:));
end


[Xmax,~] = GLgn_pts(-1*ones(1,dim),1*ones(1,dim),8);
Tineq = zeros(size(Xmax,1),length(Pf));
for r=1:1:size(Xmax,1)
    Tineq(r,:) = evaluate_MatrixOfPolys(Pf,Xmax(r,:));
end

%%
factconst = max(pn)/10;
pnfit = pn/factconst;
logpnfit = log(pnfit);

% [~,ind]=sort(pnfit(:),1,'descend');
% logppnfit=log(pnfit(ind));
% AAn=A(ind,:);
% logppnfit=logppnfit(:);



lamdim=length(lam);
K = (min(logpnfit)-5)*ones(size(Dineq,1),1);
KK=K;
DD=Dineq;
Atop=A(1:15,:);
logpntop = logpnfit(1:15);
lenconstr = length(logpnfit);

% %working good
%     minimize( 10*norm(lam2,1)+50*norm(t,2)+150*norm(t2,2))
CC=[0.5];
LAMS=zeros(lamdim,length(CC));
costs = zeros(1,length(CC));
for ci = 1:length(CC)
    cvx_begin
    variables t2(15) t(lenconstr) lam2(lamdim)
    minimize( CC(ci)*norm(lam2,1)+50*norm(t,2)+70*norm(t2,2))
    subject to
    DD*lam2 <= KK
    A*lam2==logpnfit+t
    Atop*lam2==logpntop+t2
    Tineq*lam2<=max(logpnfit)
    cvx_end
    LAMS(:,ci)=lam2;
    costs(ci) = norm(t,2);
end
[~,bind] = min(costs);
lam2 = LAMS(:,bind);
lam2(1) = lam2(1)+log(factconst);
lamsol = lam2;

% keyboard

%% normalizing and constructing normalized pdf


mxentpoly_norm=zeros(1,dim+1);
for i=1:length(lamsol)
    mxentpoly_norm=add_sub_polyND(mxentpoly_norm, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
end
mxentpoly_norm=simplify_polyND(mxentpoly_norm);

pdfnorm.func=@(x)exp(evaluate_polyND(mxentpoly_norm,x));
pdfnorm.polyeval=@(x)evaluate_polyND(mxentpoly_norm,x);
pdfnorm.poly=mxentpoly_norm;
pdfnorm.type = 'true-0I-hypercube-11';


pdfnorm = normalize_exp_pdf(pdfnorm,Xn,mquad,Pquad,'GMM_MC');


pdfnorm.trueX2normX = @(x)affineTransform(x,Atransf*Psqrt_inv,-Atransf*Psqrt_inv*mquad(:)+mulin(:));
pdfnorm.normX2trueX = @(xn)affineTransform(xn,Psqrt*inv(Atransf),mquad(:)-Psqrt*inv(Atransf)*mulin(:));
pdfnorm.normprob2trueprob = @(p)p/det(Psqrt*inv(Atransf));
pdfnorm.trueprob2normprob = @(p)p*det(Psqrt*inv(Atransf));
pdfnorm.minx = mn;
pdfnorm.maxx = mx;
pdfnorm.Pquad = Pquad;
pdfnorm.mquad = mquad;
% pdftransF.Xnorm0I2Xnorm11 = @(xn)Xnorm0I2Xnorm11(xn,minx,maxx,mquad,Pquad);
% pdftransF.Xnorm112Xnorm0I = @(xn)Xnorm112Xnorm0I(xn,minx,maxx,mquad,Pquad);



%% [pltiing and testing
% keyboard

Xt=[];
[IDX,C] = kmeans(Xn, 3);
for i=1:size(C,1)
    if length(logpn(IDX==i))>dim*2
        [m,pR]=MeanCov(Xn(IDX==i,:),pn(IDX==i)/sum(pn(IDX==i)));
        if all(eig(pR)>0)
            Xt=[Xt;mvnrnd(m,4^2*pR,200)];
        end
    end
end
NMC=size(Xt,1);
% NMC = 3000;
% Xt = [mvnrnd(zeros(dim,1),0.001^2*eye(dim),NMC/3);mvnrnd(zeros(dim,1),0.1^2*eye(dim),NMC/3);mvnrnd(zeros(dim,1),1.5^2*eye(dim),NMC/3)];
% At=zeros(NMC,length(Pf));
% tic
% for r=1:1:NMC
%    At(r,:) = evaluate_MatrixOfPolys(Pf,Xt(r,:));
% end


% lgpt=At*lamsol;
try
    lgpt=pdfnorm.polyeval(Xt);
catch
    keyboard
end
lgbnd = pdfnorm.polyeval(Xbnd);
lgpnest=pdfnorm.polyeval(Xn);
% lgpt(lgpt>max(logpn))=-10;
figure(33)
plot3(Xn(:,1),Xn(:,2),log(pn),'ro',Xn(:,1),Xn(:,2),lgpnest,'b+',Xt(:,1),Xt(:,2),lgpt,'gs')
title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
figure(34)
plot3(Xn(:,1),Xn(:,2),pn,'ro',Xn(:,1),Xn(:,2),exp(lgpnest),'b+',Xt(:,1),Xt(:,2),exp(lgpt),'gs',Xbnd(:,1),Xbnd(:,2),exp(lgbnd),'k*')
title(['time step = ',num2str(Tk),' cond = ',num2str(cond(Pquad))])
% plot_nsigellip(,1,'r',2);

%% marginal 0I 2D plots
% keyboard


% Xtrasf=pdfnorm.normX2trueX(Xn) ;
%
% figure
% plot(X(:,1),X(:,2),'ro',Xtrasf(:,1),Xtrasf(:,2),'b+')

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
% disp('Done marg')

%%

% mxentpoly=linear_transform_poly(mxentpoly_norm,Psqrt_inv,-Psqrt_inv*m(:));
%
%
% c=1/det(Psqrt);
% cexp = -log(det(Psqrt));
% c0 = get_coeff_NDpoly(mxentpoly,zeros(1,dim));
% mxentpoly = update_or_insert_coeff_NDpoly(mxentpoly,zeros(1,dim),c0+cexp);
%
% pdf.func=@(x)exp(evaluate_polyND(mxentpoly,x));
% pdf.poly=mxentpoly;


disp('Debug Stats-----')
stats.k=Tk;
stats.detPsqrt=det(Psqrt);
stats.min_pn=min(pn);
stats.max_pn=max(pn);
stats.cond=cond(Pquad);
stats


% keyboard

