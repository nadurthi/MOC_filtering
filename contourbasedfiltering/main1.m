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
time.dt=1*60*60 *constants.trueT2normT;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

%% models

model.f=@(dt,tk,xk)processmodel_2body_2D(dt,1,tk,xk);
model.fn=4;

model.h=@(x)radmodel_2D(x,1);
model.hn=1;
% model.R=diag([(0.1/constants.Re)^2,(0.5*pi/180)^2]);
model.R=diag([(0.1/constants.Re)^2]);
model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);


%% generate truth

x0=[7000,5000,1.0,3.4]';
[ r, v, Ehat ] = FnG(0, time.dt, [x0(1:2);0], [x0(3:4);0], 1);
[ r1, v1, Ehat ] = FnG(0, time.dt, r, -v, 1);

P0=diag([1^2,1^2,0.0001^2,0.0001^2]);
x0(1:2)=x0(1:2)*constants.trueX2normX;
x0(3:4)=x0(3:4)*constants.trueV2normV;
P0(1:2,1:2)=P0(1:2,1:2)*constants.trueX2normX^2;
P0(3:4,3:4)=P0(3:4,3:4)*constants.trueV2normV^2;

Xtruth = zeros(time.Ntsteps,model.fn);
Xtruth(1,:)=x0;
for k=2:time.Ntsteps
    Xtruth(k,:)=model.f(time.dt,time.Tvec(k),Xtruth(k-1,:));
end
sum(sqrt(sum(Xtruth(:,1:2).^2,2))<1)

plot(Xtruth(:,1),Xtruth(:,2),'ro')
axis equal


% plotting the propagatin of MC
Nmc=2000;
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
%    plot(XMC(:,1,k),XMC(:,2,k),'ro',Xtruth(:,1),Xtruth(:,2),'bo')
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

Npf = 5000; %paricle filter points


% generate points on contours for characterisitc solutions

Nchpol = 50;  % points used by characteristic points and polynomials

dirnmat = mvnrnd(zeros(1,model.fn),eye(model.fn),10000);
dirnmat = dirnmat./sqrt(sum(dirnmat.^2,2));
[idx,C] = kmeans(dirnmat,Nchpol);
C = C./sqrt(sum(C.^2,2));

Sphere4Dpoints = sphere4Dm(6);

Dirmats{1}=Sphere4Dpoints;

dirnmat = mvnrnd(zeros(1,model.fn),eye(model.fn),10000);
dirnmat = dirnmat./sqrt(sum(dirnmat.^2,2));
[idx,C] = kmeans(dirnmat,Nchpol);
C = C./sqrt(sum(C.^2,2));

Sphere4Dpoints = sphere4Dm(5);
Dirmats{2}=3*Sphere4Dpoints;

dirnmat = mvnrnd(zeros(1,model.fn),eye(model.fn),10000);
dirnmat = dirnmat./sqrt(sum(dirnmat.^2,2));
[idx,C] = kmeans(dirnmat,Nchpol);
C = C./sqrt(sum(C.^2,2));

Sphere4Dpoints = sphere4Dm(4);
Dirmats{3}=6*Sphere4Dpoints;

X=[Dirmats{1};Dirmats{2};Dirmats{3}]; %1sigma, 3 sigma and 6sigma
[X,w] = GH_points(zeros(4,1),0.5^2*eye(4),5);

% [X,w] = mvnrnd(zeros(4,1),0.5*eye(4),500);

A=sqrt(Pf0);
for i=1:size(X,1)
    X(i,:) = A*X(i,:)'+xf0(:);
end
probs = mvnpdf(X,xf0(:)',Pf0);

figure
plot3(X(:,1),X(:,2),X(:,3),'r+')

figure
plot(X(:,1),X(:,2),'r+')


Xinitial = X;
probsinitial = probs;

model.quadfunc=@(x,P)UT_sigmapoints(x,P,2);

[Xquad_initial,wquad_initial]=model.quadfunc(xf0(:),Pf0);
probs_quad = mvnpdf(Xquad_initial,xf0(:)',Pf0);



%% run filter
% close all

X = Xinitial;
probs = probsinitial;

Xquad=Xquad_initial;
wquad=wquad_initial;

meas_freq_steps = 5000000;

histXprior=cell(length(time.Ntsteps),5);
histXpost=cell(length(time.Ntsteps),5);

histXprior{1,1} = X;
histXprior{1,2} = probs;

histXpost{1,1} = X;
histXpost{1,2} = probs;
teststeps = [24,25];

for k=2:time.Ntsteps
    k
    Xmctest = zeros(size(XMC,1),model.fn);
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
    disp([' k = ',num2str(k)])
    
    [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(k),model);
    [Xquad,wquad]=propagate_character(Xquad,wquad,time.dt,time.Tvec(k),model);
    
    [mX,PX]=MeanCov(Xquad,wquad);
    disp(['cond = ',num2str(cond(PX))])
%         if any(k==teststeps)
    fullnormpdf=get_interp_pdf_0I(X,probs,mX,PX,4,k,Xmctest);%Xtruth(k,:)
%         end
    %     [fullpdf,pdftransF]=get_interp_pdf_hypercube11(X,probs,mX,PX,4,k,Xmctest);
    
    
    %
    %     figure(11)
    %     plot3(X(:,1),X(:,2),probs,'ro',X(:,1),X(:,2),pX,'b+',Xtestmc(:,1),Xtestmc(:,2),pXtest,'gs')
    %     title(['k = ',num2str(k)])
    
    
    histXprior{k,1}=X;
    histXprior{k,2}=probs;
    histXprior{k,3}=fullnormpdf;
    
    histXprior{k,4}=Xquad;
    histXprior{k,5}=wquad;
    
    pause(1)
    
    
    
    
    % do measurement update
    if k>=2
        if rem(k,meas_freq_steps)==0
            disp("doing meas update")
            zk = model.h(Xtruth(k,:)')+sqrtm(model.R)*randn(model.hn,1);
            zk
            
            % %%%%%%%%%%%%%%%%% MEAS UPDATE %%%%%%%%%%%%%%%%%%%%%%
            [X,probs,Xquad,wquad,fullnormpdf]=MeasUpdt_character(fullnormpdf,X,probs,Xquad,wquad,4,k,zk,model,Xtruth(k,:));
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            histXpost{k,1}=X;
            histXpost{k,2}=probs;
            histXpost{k,3}=fullnormpdf;

            histXpost{k,4}=Xquad;
            histXpost{k,5}=wquad;
    
            
        end
    end
    

    
    y=fullnormpdf.trueX2normX(X);
    py=fullnormpdf.func(y);
    probsXest=fullnormpdf.normprob2trueprob(py);
    
    figure(49)
    plot3(X(:,1),X(:,2),probs,'ro',X(:,1),X(:,2),probsXest,'b+')
    
        if k==23 %any(k==teststeps)
    
    keyboard
        end
    
    
end

% save('sim1.mat')