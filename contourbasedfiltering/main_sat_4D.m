%%
% Satellite problem

clc
close all
clear

format longg

digits(50)
%% constants
constants.radii=[6378.137,6378.137,6378.137];

constants.mu      = 3.986004418e5;     % Gravitational Const
constants.Re      = constants.radii(1);          % Earth radius (km)

constants.g0   = 9.8065;            % Sea-level acceleration, (m/s^2)

% Canonical Units
constants.muCan   = 1;
constants.RU      = constants.Re;
constants.TU      = sqrt(constants.RU^3 / constants.mu);
constants.VU      = constants.RU/constants.TU;

constants.trueA2normA=(constants.TU^2/constants.RU);
constants.normA2trueA=(constants.RU/constants.TU^2);

constants.trueV2normV=(constants.TU/constants.RU);
constants.normV2trueV=(constants.RU/constants.TU);

constants.trueX2normX=(1/constants.RU);
constants.normX2trueX=(constants.RU);

constants.trueT2normT=(1/constants.TU);
constants.normT2trueT=(constants.TU);

%% time


time.t0=0 *constants.trueT2normT;
time.tf=48*60*60 *constants.trueT2normT;
time.dt=6*60*60 *constants.trueT2normT;
time.dtplot=0.1*60*60 *constants.trueT2normT;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

time.TvecPlots=time.t0: time.dtplot :time.tf;

%% models

model.f=@(dt,tk,xk)processmodel_2body_2D(dt,1,tk,xk);
model.fn=4;
model.Q = diag([0.00001^2,0.00001^2,0.0000001^2,0.0000001^2]);
model.Q(1:2,1:2)=model.Q(1:2,1:2)*constants.trueX2normX^2;
model.Q(3:4,3:4)=model.Q(3:4,3:4)*constants.trueV2normV^2;
model.fstates = {'x','y','vx','vy'};


model.h=@(xnorm)radmodel_2D(xnorm,constants);
model.hn=1;
model.hvec=@(x)sensmodel2vec(x,model.h,model.hn);
% model.R=diag([(2*pi/180)^2]);
model.R=diag([(0.5*pi/180)^2]);
model.hstates = {'azi'};
% model.R=diag([(0.1/constants.Re)^2]);
model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);



%% generate truth

x0=[10000,10,-0.1,7.6]';
OE = cart2orbelem(sat4D26D(x0),constants.mu)
% [ r, v, Ehat ] = FnG(0, time.dt, x0(1:3), x0(4:6), 1);
% [ r1, v1, Ehat ] = FnG(0, time.dt, r, -v, 1);

P0=diag([1^2,1^2,0.001^2,0.001^2]);

x0(1:2)=x0(1:2)*constants.trueX2normX;
x0(3:4)=x0(3:4)*constants.trueV2normV;
P0(1:2,1:2)=P0(1:2,1:2)*constants.trueX2normX^2;
P0(3:4,3:4)=P0(3:4,3:4)*constants.trueV2normV^2;

Xtruth = zeros(time.Ntsteps,model.fn);
Xtruth(1,:)=x0;
for k=2:time.Ntsteps
    Xtruth(k,:)=model.f(time.dt,time.Tvec(k-1),Xtruth(k-1,:));
end
Xplot = zeros(time.Ntsteps,model.fn);
Xplot(1,:)=x0;
for k=2:length(time.TvecPlots)
    Xplot(k,:)=model.f(time.dtplot,time.TvecPlots(k-1),Xplot(k-1,:));
end

figure
plot(Xtruth(:,1),Xtruth(:,2),'ro')
hold on
plot(Xplot(:,1),Xplot(:,2),'b')
axis equal
axis square

% plotting the propagatin of MC
Nmc=5000;
XMC=zeros(Nmc,model.fn,time.Ntsteps);
XMC(:,:,1)=mvnrnd(x0',P0,Nmc);
for i=1:Nmc
    i
    for k=2:time.Ntsteps
        XMC(i,:,k)=model.f(time.dt,time.Tvec(k-1),XMC(i,:,k-1));
    end
end
pMC = mvnpdf(XMC(:,:,1),x0',P0);

%%
close all
if false
    for k=1:time.Ntsteps
        figure(1)
        clf
        subplot(1,2,1)
        plot3(XMC(:,1,k),XMC(:,2,k),log(pMC),'ro')
        hold on
        plot3(Xplot(:,1),Xplot(:,2),Xplot(:,3),'b-')
        title(['k = ',num2str(k)])
        axis equal
        axis square
        
        subplot(1,2,2)
        XXX = zeros(Nmc,model.fn);
        XXX = XMC(:,:,k);
        [m,P] = MeanCov(XXX,ones(Nmc,1)/Nmc);
        A=sqrtm(inv(P));
        for i=1:Nmc
            XXX(i,:) = A*(XXX(i,:)'-m(:));
        end
        plot3(XXX(:,1),XXX(:,2),log(pMC),'ro')
        title(['k = ',num2str(k)])
        axis equal
        axis square
        view([-66,90])
        keyboard
    end
end
%% comparing with UKF and particle filter
xf0=mvnrnd(x0(:)',P0);
% xf0 = x0;
Pf0 = P0;

xfquad = xf0;
Pfquad = P0;

Npf = 5000; %paricle filter points


GMMinitial=getInitialGMM(xf0(:),0.9^2*Pf0,7);
model.gmmmethod='ut';

% generate points on contours for characterisitc solutions

model.pointGenerator = @(mx,Px)GH_points(mx,Px,9);
% model.pointGenerator = @(mx,Px)mvnrnd(mx(:)',Px,20000);

model.pointscalefactor = 0.4;
[X,w] = model.pointGenerator(zeros(model.fn,1),model.pointscalefactor^2*eye(model.fn));

% [X1,w] = GH_points(zeros(model.fn,1),0.5^2*eye(model.fn),5);
% % [X2,w] = GH_points(zeros(model.fn,1),2^2*eye(model.fn),4);
% X2 = 6*sphere6Dm(6); %3906
% % [X,w] = mvnrnd(zeros(4,1),0.5*eye(4),500);
% X=[X1;X2];

A=sqrt(Pf0);
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
probs_quad = mvnpdf(Xquad_initial,xf0(:)',Pf0);



%% run filter
close all
clc

X = Xinitial;
probs = probsinitial;

Xquad=Xquad_initial;
wquad=wquad_initial;

xfquad = xf0;
Pfquad = P0;

GMM = GMMinitial;


meas_freq_steps = 1;

histXprior=cell(time.Ntsteps,5);
histXpost=cell(time.Ntsteps,5);

histXprior{1,1} = X;
histXprior{1,2} = probs;

histXpost{1,1} = X;
histXpost{1,2} = probs;
teststeps = [33];

plotfolder='simulations/SAT4meas_test2';
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

[mX,PX]=MeanCov(X,probs/sum(probs));
EstMOCfilter_mu(1,:) = mX;
EstMOCfilter_P(1,:) = reshape(PX,1,model.fn^2);

EstQuadfilter_mu(1,:) = xfquad;
EstQuadfilter_P(1,:) = reshape(Pfquad,1,model.fn^2);


for k=2:time.Ntsteps
%     close all
    
    k
    Xmctest = zeros(size(XMC,1),model.fn);
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
    Xmctest=[];
    disp([' k = ',num2str(k)])
    
    [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(k),model);
    [Xquad,wquad]=propagate_character(Xquad,wquad,time.dt,time.Tvec(k),model);
    [xfquad,Pfquad]=QuadProp(xfquad,Pfquad,time.dt,time.Tvec(k),model,'ut');

%     GMM=prior_prop_GMM(GMM,time.dt,time.Tvec(k),model);
    
    
    [mX,PX]=MeanCov(X,probs/sum(probs));
    

    disp(['cond = ',num2str(cond(PX))])

    
%         if any(k==teststeps)
%     keyboard
%     fullnormpdf=get_interp_pdf_0I_boostmixGaussian(X,probs,mX,PX,3,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)
%     fullnormpdf=get_interp_pdf_0I_2D(X,probs,mX,PX,4,k,[],Xtruth(k,:),plotsconf); %Xtruth(k,:)
    fullnormpdf=get_interp_pdf_0I_4Dsat(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)
%     fullnormpdf=get_interp_pdf_0I_2D(X,probs,mX,PX,4,k,[],Xtruth(k,:),plotsconf); %Xtruth(k,:)
    priorpdfnorm=fullnormpdf;
    
    figure(1)
    plotpdfs_prior_4D([1,2],k,fullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
    
%     keyboard
    
    histXprior{k,1}=X;
    histXprior{k,2}=probs;
    histXprior{k,3}=fullnormpdf;
    
    histXprior{k,4}=xfquad;
    histXprior{k,5}=Pfquad;
    
%    keyboard
    pause(1)
    
    
    
    
    % do measurement update
    if k>=2
        if rem(k,meas_freq_steps)==0
            disp("doing meas update")
            zk = model.h(Xtruth(k,:)')+sqrtm(model.R)*randn(model.hn,1);
            zk
            
            % %%%%%%%%%%%%%%%%% MEAS UPDATE %%%%%%%%%%%%%%%%%%%%%%
%             [X,probs,fullnormpdf] =       MeasUpdt_character_modf(X,probs,3,k,zk,Xtruth(k,:),model,Xmctest,11);
%             keyboard
            [X,probs,fullnormpdf] = MeasUpdt_character_modf_4Dsat(X,probs,priorpdfnorm,4,k,zk,Xtruth(k,:),model,Xmctest,5);
            
%             [X,probs,fullnormpdf]=MeasUpdt_character_modf(fullnormpdf,X,probs,4,k,zk,Xtruth,model,Xmctest);
            [xfquad,Pfquad]=QuadMeasUpdt(xfquad,Pfquad,zk,time.dt,time.Tvec(k),model,'ut');
            
            GMM=post_GMM(GMM,zk,time.dt,time.Tvec(k),model);
            
            
            plotpdfs_post_4D([1,2],k,priorpdfnorm,fullnormpdf,model,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),zk,savePostProps)
    
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            histXpost{k,1}=X;
            histXpost{k,2}=probs;
            histXpost{k,3}=fullnormpdf;

            histXpost{k,4}=xfquad;
            histXpost{k,5}=Pfquad;
    
            
        end
    end

%     keyboard

%     pause(1)
    
    
    [mestX,PestX]=MeanCov(X,probs/sum(probs));

    
    EstMOCfilter_mu(k,:) = mestX;
    EstMOCfilter_P(k,:) = reshape(PestX,1,model.fn^2);

    EstQuadfilter_mu(k,:) = xfquad;
    EstQuadfilter_P(k,:) = reshape(Pfquad,1,model.fn^2);

end

save([plotfolder,'/sim1.mat'])