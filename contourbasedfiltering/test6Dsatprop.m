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
time.dt=16*60*60 *constants.trueT2normT;
time.dtplot=0.1*60*60 *constants.trueT2normT;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

time.TvecPlots=time.t0: time.dtplot :time.tf;

%% models

model.f=@(dt,tk,xk)processmodel_2body(dt,1,tk,xk);
model.fback=@(dt,tk,xk)processmodel_2body_backward(dt,1,tk,xk);

% [a,e,i,om,Om,M]
model.fn=6;
model.Q = diag([0.00001^2,0.00001^2,0.00001^2,0.0000001^2,0.0000001^2,0.0000001^2]);
model.Q(1:3,1:3)=model.Q(1:3,1:3)*constants.trueX2normX^2;
model.Q(4:6,4:6)=model.Q(4:6,4:6)*constants.trueV2normV^2;
% model.fstates = {'a','e','i','om','Om','M'};
model.fstates = {'x','y','z','vx','vy','vz'};

model.h=@(x)radmodel(xk,1);
model.hn=2;
model.hvec=@(x)sensmodel2vec(x,model.h,model.hn);
model.R=diag([(0.1*pi/180)^2,(0.1*pi/180)^2]);
model.hstates = {'azi','el'};
% model.R=diag([(0.1/constants.Re)^2]);
model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);



%% generate truth

x0=[10000,10,3000,-0.1,7.4,2]';
OE = cart2orbelem(x0,constants.mu)
% [ r, v, Ehat ] = FnG(0, time.dt, x0(1:3), x0(4:6), 1);
% [ r1, v1, Ehat ] = FnG(0, time.dt, r, -v, 1);

P0=diag([1^2,1^2,1^2,0.001^2,0.001^2,0.001^2]);

x0(1:3)=x0(1:3)*constants.trueX2normX;
x0(4:6)=x0(4:6)*constants.trueV2normV;
P0(1:3,1:3)=P0(1:3,1:3)*constants.trueX2normX^2;
P0(4:6,4:6)=P0(4:6,4:6)*constants.trueV2normV^2;



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



xf0 = x0;
Pf0 = P0;

xfquad = xf0;
Pfquad = P0;

Npf = 5000; %paricle filter points


GMMinitial=getInitialGMM(xf0(:),0.9^2*Pf0,4);
model.gmmmethod='ut';

% generate points on contours for characterisitc solutions

model.pointGenerator = @(mx,Px)GH_points(mx,Px,3);
% model.pointGenerator = @(mx,Px)mvnrnd(mx(:)',Px,20000);

model.pointscalefactor = 0.9;
[Xut,wut] = model.pointGenerator(zeros(model.fn,1),model.pointscalefactor^2*eye(model.fn));

model.pointscalefactor = 0.9;
[X,w] = model.pointGenerator(zeros(model.fn,1),model.pointscalefactor^2*eye(model.fn));


pf.no_particles =100;
pf.no_bins = 100;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;         
pf.regularizepf = 1;


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


Npt=size(Xinitial,1);
Xpttraj=zeros(Npt,model.fn,time.Ntsteps);
Xpttraj(:,:,1)=Xinitial;
for i=1:Npt
    i
    for k=2:time.Ntsteps
        Xpttraj(i,:,k)=model.f(time.dt,time.Tvec(k-1),Xpttraj(i,:,k-1));
    end
end
p=4;
Xr = zeros(Npt,model.fn);
for ii=1:Npt
    Xr(ii,:) = Xpttraj(ii,:,p);
end

figure
plot3(Xr(:,1),Xr(:,2),Xr(:,3),'bo')





%% run filter
close all
clc

X = Xinitial;
probs = probsinitial;

Xquad=Xquad_initial;
wquad=wquad_initial;

xfquad = xf0;
Pfquad = P0;

dsX = DataSet(X,probs,'TrueState');
dsX.AddMeanCov_to_OI_Trasform(xf0(:),Pf0);

priorfullnormpdf.func=@(x)mvnpdf(x,zeros(model.fn,1)',eye(model.fn));
priorfullnormpdf.transForms = dsX.GetTrasnformers();

priorfullnormpdf.info = 'true-0I-hypercube-11';
priorfullnormpdf.pdftype = 'ExpPdf';

priorfullnormpdf.GMMHull = NaN;
priorfullnormpdf.LB = NaN;
priorfullnormpdf.UB = NaN;
priorfullnormpdf.RigTree = NaN;

% initial condition for the filter
% X_pf=repmat(mu_pf,1,pf.no_particles)+sqrtm(P_pf)*randn(model.fn,pf.no_particles);
% w_pf = ones(1, pf.no_particles) / pf.no_particles;

GMM = GMMinitial;


meas_freq_steps = 1000000;

histXprior=cell(time.Ntsteps,5);
histXpost=cell(time.Ntsteps,5);

histXprior{1,1} = X;
histXprior{1,2} = probs;

histXpost{1,1} = X;
histXpost{1,2} = probs;
teststeps = [33];

plotfolder='simulations/SAT6Dsim1_props';
mkdir(plotfolder)

savePriorProps.plotfolder=[plotfolder,'/prior'];
mkdir(savePriorProps.plotfolder)
savePriorProps.saveit=0;

savePostProps.plotfolder=[plotfolder,'/post'];
mkdir(savePostProps.plotfolder)
savePostProps.saveit=0;

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
%     Xmctest=[];
    disp([' k = ',num2str(k)])
    
    [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(k),model);
    [Xquad,wquad]=propagate_character(Xquad,wquad,time.dt,time.Tvec(k),model);
    [xfquad,Pfquad]=QuadProp(xfquad,Pfquad,time.dt,time.Tvec(k),model,'ut');

%     GMM=prior_prop_GMM(GMM,time.dt,time.Tvec(k),model);
    
%      X_pf = pftimeupdate(X_pf, model);
%      w_pf=w_pf;
 
 
    [mX,PX]=MeanCov(X,probs/sum(probs));
    

    disp(['cond = ',num2str(cond(PX))])

    priorfullnormpdf = get_interp_pdf_0I_6D_gobackget(model,X,probs,4,k,time.Tvec(k),time.dt,priorfullnormpdf,Xmctest,Xtruth(k,:));
%     fullnormpdf=get_interp_pdf_0I_6Dsat(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)

%     priorpdfnorm=fullnormpdf;
    
    
    plotpdfs_prior_6D([1,2],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
    
end