%%

clc
close all
clear

format longg

digits(50)

%% time


time.t0=0;
time.tf=10;
time.dt=0.5;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

%% models

model.f=@(dt,tk,xk)duff_prop_model(dt,xk);
model.fback=@(dt,tk,xk)duff_prop_model_back(dt,xk);

model.fn=2;
model.Q = diag([(1e-6)^2,(1e-6)^2]);
model.sQ = sqrtm(model.Q);
model.fstates = {'x1','x2','v1','v2'};

% model.h=@(x)[sqrt(x(1)^2+x(2)^2)];
model.h=@(x)duff_meas_model_th(x);
model.hn=1;
model.hvec=@(x)sensmodel2vec(x,model.h,model.hn);

% model.hstates = {'x','y'};
% model.hstates = {'r','th'};
model.hstates = {'th'};
% model.R=diag([(0.1/constants.Re)^2,(0.5*pi/180)^2]);
% model.R=diag([0.1^2,0.1^2]);
% model.R=diag([1^2,(2*pi/180)^2]);
model.R=diag([(1*pi/180)^2]);
% model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);


%% generate truth

x0=[5,5]';
P0=0.5^2*eye(2);

Xtruth = zeros(time.Ntsteps,model.fn);
Xtruth(1,:)=x0;
for k=2:time.Ntsteps
    Xtruth(k,:)=model.f(time.dt,time.Tvec(k),Xtruth(k-1,:));
end

Ntsteps=5000;
plotTvec=linspace(time.t0,time.tf,Ntsteps);
plotdt = plotTvec(2)-plotTvec(1);
Xtruplot = zeros(Ntsteps,model.fn);
Xtruplot(1,:)=x0;
for k=2:Ntsteps
    Xtruplot(k,:)=model.f(plotdt,plotTvec(k),Xtruplot(k-1,:));
end


figure
plot(Xtruth(:,1),Xtruth(:,2),'ro')
axis equal
axis square

figure
[t,x]=ode45(@duff,linspace(0,10,1000),x0);
plot(x(:,1),x(:,2),'b')
axis equal
axis square


% plotting the propagatin of MC
Nmc=5000;
XMC=zeros(Nmc,model.fn,time.Ntsteps);
XMC(:,:,1)=mvnrnd(x0',P0,Nmc);
for i=1:Nmc
    for k=2:time.Ntsteps
        XMC(i,:,k)=model.f(time.dt,time.Tvec(k),XMC(i,:,k-1));
    end
end

%%

% for k=1:time.Ntsteps
%    figure(1)
%    plot(XMC(:,1,k),XMC(:,2,k),'ro',Xtruth(:,1),Xtruth(:,2),'b')
%    title(['k = ',num2str(k)])
%    axis equal
%    axis square
% 
%    figure(2)
%    plot(XMC(:,1,k),XMC(:,2,k),'ro')
%    title(['k = ',num2str(k)])
%    axis equal
%    axis square
% 
%    pause(1)
% end

%% comparing with UKF and particle filter
xf0=mvnrnd(x0(:)',P0);
xf0 = x0;
Pf0 = P0;

xfquad = xf0;
Pfquad = P0;

% Npf = 2000; %paricle filter points


% GMMinitial=getInitialGMM(xf0(:),0.9^2*Pf0,4);
% model.gmmmethod='ut';

% generate points on contours for characterisitc solutions
model.pointscalefactor = 0.9;
model.pointGenerator = @(mx,Px)hybridpoints(mx,Px,3,0.2,5,0.6);
% model.pointGenerator = @(mx,Px)mvnrnd(mx(:)',Px,20000);


% [Xut,wut] = UT_sigmapoints(zeros(model.fn,1),eye(model.fn),2);
% 
% [X,w] = model.pointGenerator(zeros(model.fn,1),eye(model.fn));


pf.no_particles =25000;
pf.no_bins = 100;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;         
pf.regularizepf = 0;


% [X,w] = GH_points(zeros(model.fn,1),1*eye(model.fn),4);
% [X2,w] = GH_points(zeros(model.fn,1),2^2*eye(model.fn),4);
% X2 = 6*sphere6Dm(6); %3906
[X,w] = mvnrnd(zeros(model.fn,1),1*eye(model.fn),500);
% X=[X1;X2];

A=sqrtm(Pf0);
for i=1:size(X,1)
    X(i,:) = A*X(i,:)'+xf0(:);
end
probs = mvnpdf(X,xf0(:)',Pf0);


figure
plot(X(:,1),X(:,2),'r+')


Xinitial = X;
probsinitial = probs;

model.quadfunc=@(x,P)UT_sigmapoints(x,P,2);

[Xquad_initial,wquad_initial]=model.quadfunc(xf0(:),Pf0);
% probs_quad = mvnpdf(Xquad_initial,xf0(:)',Pf0);


Npt=size(Xinitial,1);
Xpttraj=zeros(Npt,model.fn,time.Ntsteps);
Xpttraj(:,:,1)=Xinitial;
for i=1:Npt
    i
    for k=2:time.Ntsteps
        Xpttraj(i,:,k)=model.f(time.dt,time.Tvec(k-1),Xpttraj(i,:,k-1));
    end
end
% p=4;
% Xr = zeros(Npt,model.fn);
% for ii=1:Npt
%     Xr(ii,:) = Xpttraj(ii,:,p);
% end
% 
% figure
% plot3(Xr(:,1),Xr(:,2),Xr(:,3),'bo')




%% run filter
close all
clc

X_pf=mvnrnd(xf0(:)',P0,pf.no_particles)';
w_pf = 1/pf.no_particles * ones(pf.no_particles,1);



xfquad = xf0;
Pfquad = P0;



% [X1,~] = GH_points(xf0(:),Pf0,4);
% probsk1 = mvnpdf(X1,xf0(:)',Pf0);
dsX = DataSet(Xinitial,probsinitial,'TrueState');
dsX.AddMeanCov_to_OI_Trasform(xf0(:),5^2*Pf0);
[ma,Pb]=MeanCov(dsX.X,dsX.p/sum(dsX.p));
% dsX.AddHyperCubeTrasform(-0.9*ones(dim,1),0.9*ones(dim,1));
transForms = dsX.GetTrasnformers();

pdfk.Xrep = Xinitial;
pdfk.prep = probsinitial;
pdfk.normeval=@(x)mvnpdf(x,ma(:)',Pb);
pdfk.transForms = dsX.GetTrasnformers();

pdfk.info = 'true-0I';
pdfk.pdftype = 'ExpPdf';

% pdfnorm.GMMHull = GMMHull;
pdfk.LB = -2*ones(1,model.fn);
pdfk.UB = 2*ones(1,model.fn);


meas_freq_steps = 1;


plotfolder='simulations/duff2D_measUPDATE_treefilter2';
mkdir(plotfolder)

savePriorProps.plotfolder=[plotfolder,'/prior'];
mkdir(savePriorProps.plotfolder)
savePriorProps.saveit=1;

savePostProps.plotfolder=[plotfolder,'/post'];
mkdir(savePostProps.plotfolder)
savePostProps.saveit=1;

EstMOCfilter_mu =   zeros(time.Ntsteps,model.fn);
EstMOCfilter_P =    zeros(time.Ntsteps,model.fn^2);

EstQuadfilter_mu =   zeros(time.Ntsteps,model.fn);
EstQuadfilter_P =    zeros(time.Ntsteps,model.fn^2);

EstPFfilter_mu =   zeros(time.Ntsteps,model.fn);
EstPFfilter_P =    zeros(time.Ntsteps,model.fn^2);

[mX,PX]=MeanCov(X,probs/sum(probs));
EstMOCfilter_mu(1,:) = mX;
EstMOCfilter_P(1,:) = reshape(PX,1,model.fn^2);

EstQuadfilter_mu(1,:) = xfquad;
EstQuadfilter_P(1,:) = reshape(Pfquad,1,model.fn^2);

[mpf,Ppf]=MeanCov(X_pf',w_pf);
EstPFfilter_mu(1,:) = mpf;
EstPFfilter_P(1,:) = reshape(Ppf,1,model.fn^2);

postfullnormpdf=pdfk;


histXprior{1,1}=postfullnormpdf;
histXprior{1,2}=xfquad;
histXprior{1,3}=Pfquad;
histXprior{1,4}=X_pf;
histXprior{1,5}=w_pf;


histXpost{1,1}=postfullnormpdf;
histXpost{1,2}=xfquad;
histXpost{1,3}=Pfquad;
histXpost{1,4}=X_pf;
histXpost{1,5}=w_pf;


postfullnormpdf 
for k=2:time.Ntsteps
%     close all
    
    k
    Xmctest = zeros(size(XMC,1),model.fn);
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
    Xmctest=[];
    
    disp([' k = ',num2str(k)])
    

    [xfquad,Pfquad]=QuadProp(xfquad,Pfquad,time.dt,time.Tvec(k),model,'ut');
    
     X_pf = pftimeupdate_custommoc(X_pf, model,time.Tvec(k),k,time.dt);
     w_pf=w_pf;
    

     
    % do measurement update
    if k>=2
        if rem(k,meas_freq_steps)==0
            disp("doing meas update")
            zk = model.h(Xtruth(k,:)')+sqrtm(model.R)*randn(model.hn,1);
            zk
            
            % %%%%%%%%%%%%%%%%% MEAS UPDATE %%%%%%%%%%%%%%%%%%%%%%
            
            postfullnormpdfkm1 = postfullnormpdf;
            [priorfullnormpdf,postfullnormpdf,Xnptk1_prior,probntk1_prior,Xnptk1_post,probntk1_post]=moc_pf(postfullnormpdf,k-1,time,zk,model,struct('genpriorpdf',1,'genpostpdf',1));
%             postfullnormpdf = priorfullnormpdf;
    
            
            xfquadprior=xfquad;
            Pfquadprior=Pfquad;
            [xfquad,Pfquad]=QuadMeasUpdt(xfquad,Pfquad,zk,time.dt,time.Tvec(k),model,'ut');
%             GMM=post_GMM(GMM,zk,time.dt,time.Tvec(k),model);
            X_pfprior=X_pf;
            w_pfprior=w_pf;
            [X_pf, w_pf] = pfmeasupdate_custommoc(X_pf, w_pf, model,time.Tvec(k),k,time.dt, pf, zk);
            

%             keyboard
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            

            
        end
    end
    
    % plot PRIOR --------------
%     plorprior_treepdf(k,model,Xtruplot,Xtruth,priorfullnormpdf,Xmctest,Xnptk1_prior,probntk1_prior,X_pf,xfquad,Pfquad,savePriorProps)
%     
    % plot POST --------------
    plotpost_treepdf(k,model,Xtruplot,Xtruth,priorfullnormpdf,postfullnormpdf,Xmctest,Xnptk1_prior,probntk1_prior,Xnptk1_post,probntk1_post,X_pf,w_pf,xfquadprior,Pfquadprior,xfquad,Pfquad,savePostProps)

    
%     keyboard
    
%     plot3(Xnptk1_prior(:,1),Xnptk1_prior(:,2),probntk1_prior,'bo')
%     hold on
%     plot3Dtreepatches(priorfullnormpdf.boxes)
    
    
    
    histXprior{k,1}=priorfullnormpdf;
    histXprior{k,2}=xfquadprior;
    histXprior{k,3}=Pfquadprior;
    histXprior{k,4}=X_pfprior;
    histXprior{k,5}=w_pfprior;

    histXpost{k,1}=postfullnormpdf;
    histXpost{k,2}=xfquad;
    histXpost{k,3}=Pfquad;
    histXpost{k,4}=X_pf;
    histXpost{k,5}=w_pf;

    [mestX,PestX]=MeanCov(postfullnormpdf.Xrep,postfullnormpdf.prep/sum(postfullnormpdf.prep));
    [mpf,Ppf]=MeanCov(X_pf',w_pf(:)/sum(w_pf));
    
    EstMOCfilter_mu(k,:) = mestX;
    EstMOCfilter_P(k,:) = reshape(PestX,1,model.fn^2);

    EstQuadfilter_mu(k,:) = xfquad;
    EstQuadfilter_P(k,:) = reshape(Pfquad,1,model.fn^2);
    
    EstPFfilter_mu(k,:) = mpf;
    EstPFfilter_P(k,:) = reshape(Ppf,1,model.fn^2);
    
end

%% Estimate plots
Err_Xpost=zeros(time.Ntsteps,model.fn);
P_Xpost=zeros(time.Ntsteps,model.fn^2);
Sig_Xpost=zeros(time.Ntsteps,model.fn);

Err_quad=zeros(time.Ntsteps,model.fn);
P_quad=zeros(time.Ntsteps,model.fn^2);
Sig_quad=zeros(time.Ntsteps,model.fn);

for k=2:time.Ntsteps
    X=histXpost{k,1};
    probs=histXpost{k,2};
    
    xfquad = histXpost{k,4};
    Pfquad = histXpost{k,5};
    
    [mestX,PestX]=MeanCov(X,probs/sum(probs));
    
    Err_Xpost(k,:) = Xtruth(k,:)-mestX(:)';
    P_Xpost(k,:) = reshape(PestX,1,model.fn^2);
    Sig_Xpost(k,:) = diag(sqrtm(PestX));
    
    Err_quad(k,:) = Xtruth(k,:)-xfquad(:)';
    P_quad(k,:) = reshape(Pfquad,1,model.fn^2); 
    Sig_quad(k,:) = diag(sqrtm(Pfquad));
    
end

fff=plotfolder;
fff='simulations/duffsim_thmeas_polyGHGMMresample_good';
ppp=1;
figure
plot(time.Tvec(2:end),Err_Xpost(2:end,ppp),'b',time.Tvec(2:end),Err_Xpost(2:end,ppp)+ 3*Sig_Xpost(2:end,ppp),'r--',time.Tvec(2:end),Err_Xpost(2:end,ppp)- 3*Sig_Xpost(2:end,ppp),'r--','linewidth',2)
xlabel('time')
ylabel('x_1 - error')
grid
plot_prop_paper_2D
saveas(gcf,[fff,'/','reg_x1err'],'png')
saveas(gcf,[fff,'/','reg_x1err'],'fig')

ppp=1;
figure
plot(time.Tvec(2:end),Err_quad(2:end,ppp),'b',time.Tvec(2:end),Err_quad(2:end,ppp)+ 3*Sig_quad(2:end,ppp),'r--',time.Tvec(2:end),Err_quad(2:end,ppp)- 3*Sig_quad(2:end,ppp),'r--','linewidth',2)
xlabel('time')
ylabel('x_1 - error')
grid
plot_prop_paper_2D
saveas(gcf,[fff,'/','quad_x1err'],'png')
saveas(gcf,[fff,'/','quad_x1err'],'fig')

ppp=2;
figure
plot(time.Tvec(2:end),Err_Xpost(2:end,ppp),'b',time.Tvec(2:end),Err_Xpost(2:end,ppp)+ 3*Sig_Xpost(2:end,ppp),'r--',time.Tvec(2:end),Err_Xpost(2:end,ppp)- 3*Sig_Xpost(2:end,ppp),'r--','linewidth',2)
xlabel('time')
ylabel('x_2 - error')
grid
plot_prop_paper_2D
saveas(gcf,[fff,'/','reg_x2err'],'png')
saveas(gcf,[fff,'/','reg_x2err'],'fig')

ppp=2;
figure
plot(time.Tvec(2:end),Err_quad(2:end,ppp),'b',time.Tvec(2:end),Err_quad(2:end,ppp)+ 3*Sig_quad(2:end,ppp),'r--',time.Tvec(2:end),Err_quad(2:end,ppp)- 3*Sig_quad(2:end,ppp),'r--','linewidth',2)
xlabel('time')
ylabel('x_2 - error')
grid
plot_prop_paper_2D
saveas(gcf,[fff,'/','quad_x2err'],'png')
saveas(gcf,[fff,'/','quad_x2err'],'fig')




%%
save([plotfolder,'/sim1.mat'])