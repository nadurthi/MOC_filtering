%%
% Satellite problem

clc
close all
clear

format longg

digits(50)

%% time


time.t0=0;
time.tf=10;
time.dt=0.2;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

%% models

model.f=@(dt,tk,xk)duff_prop_model(dt,xk);
model.fn=2;
model.Q = diag([(1e-4)^2,(1e-4)^2]);


model.h=@(x)[sqrt(x(1)^2+x(2)^2)];
model.hn=1;
% model.R=diag([(0.1/constants.Re)^2,(0.5*pi/180)^2]);
model.R=diag([0.5^2]);
model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);


%% generate truth

x0=[5,5]';
P0=0.1^2*eye(2);

Xtruth = zeros(time.Ntsteps,model.fn);
Xtruth(1,:)=x0;
for k=2:time.Ntsteps
    Xtruth(k,:)=model.f(time.dt,time.Tvec(k),Xtruth(k-1,:));
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
Nmc=1000;
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
% xf0=mvnrnd(x0(:)',P0);
xf0 = x0;
Pf0 = P0;

xfquad = xf0;
Pfquad = P0;

Npf = 5000; %paricle filter points


% generate points on contours for characterisitc solutions

model.pointGenerator = @(mx,Px)GH_points(mx,Px,11);
[X,w] = model.pointGenerator(zeros(model.fn,1),0.5^2*eye(model.fn));


% [X,w] = mvnrnd(zeros(4,1),0.5*eye(4),500);

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

X = Xinitial;
probs = probsinitial;

Xquad=Xquad_initial;
wquad=wquad_initial;

xfquad = xf0;
Pfquad = P0;

meas_freq_steps = 1;

histXprior=cell(time.Ntsteps,5);
histXpost=cell(time.Ntsteps,5);

histXprior{1,1} = X;
histXprior{1,2} = probs;

histXpost{1,1} = X;
histXpost{1,2} = probs;
teststeps = [33];

plotfolder='duffsim1_prop';
mkdir(plotfolder)

savePriorProps.plotfolder=plotfolder;
savePriorProps.saveit=0;

savePostProps.plotfolder=plotfolder;
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
    Xmctest=[];
    disp([' k = ',num2str(k)])
    
    [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(k),model);
    [Xquad,wquad]=propagate_character(Xquad,wquad,time.dt,time.Tvec(k),model);
    [xfquad,Pfquad]=QuadProp(xfquad,Pfquad,time.dt,time.Tvec(k),model,'ut');

    
    
    [mX,PX]=MeanCov(X,probs/sum(probs));
    

    disp(['cond = ',num2str(cond(PX))])

    
%         if any(k==teststeps)
    
    fullnormpdf=get_interp_pdf_0I_boostmixGaussian(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)
%     fullnormpdf=get_interp_pdf_0I_2D(X,probs,mX,PX,4,k,[],Xtruth(k,:),plotsconf); %Xtruth(k,:)
    
    
    plotpdfs_prior_2D(k,fullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),savePriorProps)
    
    histXprior{k,1}=X;
    histXprior{k,2}=probs;
    histXprior{k,3}=fullnormpdf;
    
    histXprior{k,4}=xfquad;
    histXprior{k,5}=Pfquad;
    
   
    pause(1)
    
    
    
    
    % do measurement update
    if k>=2
        if rem(k,meas_freq_steps)==0
            disp("doing meas update")
            zk = model.h(Xtruth(k,:)')+sqrtm(model.R)*randn(model.hn,1);
            zk
            
            % %%%%%%%%%%%%%%%%% MEAS UPDATE %%%%%%%%%%%%%%%%%%%%%%
            [X,probs,fullnormpdf]=MeasUpdt_character_modf(X,probs,4,k,zk,Xtruth(k,:),model,Xmctest,11);
%             [X,probs,fullnormpdf]=MeasUpdt_character_modf(fullnormpdf,X,probs,4,k,zk,Xtruth,model,Xmctest);
            [xfquad,Pfquad]=QuadMeasUpdt(xfquad,Pfquad,zk,time.dt,time.Tvec(k),model,'ut');
            
            plotpdfs_post_2D(k,fullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),savePostProps)
    
            
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
    [mestquad,PestXquad]=MeanCov(xfquad,Pfquad);
    
    EstMOCfilter_mu(1,:) = mestX;
    EstMOCfilter_P(1,:) = reshape(PestX,1,model.fn^2);

    EstQuadfilter_mu(1,:) = mestquad;
    EstQuadfilter_P(1,:) = reshape(PestXquad,1,model.fn^2);

end

% save('sim1.mat')