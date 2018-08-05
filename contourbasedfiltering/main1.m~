%%
% Satellite problem

clc
close all
clear

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
time.dt=5*60 *constants.trueT2normT;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

%% models

model.f=@(dt,tk,xk)processmodel_2body_2D(dt,1,tk,xk);
model.fn=4;

model.h=@(x)radmodel_2D(x);
model.hn=2;
model.R=diag([0.1^2,(5*pi/180)^2]);
model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);


%% generate truth

x0=[7000,5000,2.0,4.4]';
[ r, v, Ehat ] = FnG(0, time.dt, [x0(1:2);0], [x0(3:4);0], 1);
[ r1, v1, Ehat ] = FnG(0, time.dt, r, -v, 1);

P0=diag([1^2,1^2,0.01^2,0.01^2]);
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
Nmc=100;
XMC=zeros(Nmc,model.fn,time.Ntsteps);
XMC(:,:,1)=mvnrnd(x0',P0,Nmc);
for i=1:Nmc
    for k=2:time.Ntsteps
        XMC(i,:,k)=model.f(time.dt,time.Tvec(k),XMC(i,:,k-1));
    end
end

% figure
% for k=1:time.Ntsteps
% 
%    plot(XMC(:,1,k),XMC(:,2,k),'ro')
%    title(['k = ',num2str(k)])
%    axis equal
%    axis square
%    pause(1)
% end

%% comparing with UKF and particle filter

xf0 = mvnrnd(x0,P0);
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
% [X,w] = GH_points(zeros(4,1),eye(4),5);

[X,w] = mvnrnd(zeros(4,1),eye(4),500);

A=sqrt(Pf0);
for i=1:size(X,1)
    X(i,:) = A*X(i,:)'+xf0';
end
probs = mvnpdf(X,xf0,Pf0);

figure
plot3(X(:,1),X(:,2),X(:,3),'r+')

figure
plot(X(:,1),X(:,2),'r+')


Xinitial = X;
probsinitial = probs;

[Xquad,wquad]=UT_sigmapoints(xf0,Pf0,2);
probs_quad = mvnpdf(Xquad,xf0,Pf0);



%% run filter
close all

X = Xinitial;
probs = probsinitial;

meas_freq_steps = 100000;
histXprior=cell(length(time.Ntsteps),2);
histXpost=cell(length(time.Ntsteps),2);

histXprior{1,1} = X;
histXprior{1,2} = probs;

histXpost{1,1} = X;
histXpost{1,2} = probs;

for k=2:time.Ntsteps
    k
    disp([' k = ',num2str(k)])
    
    [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(k),model);
    [Xquad,wquad]=propagate_character(Xquad,wquad,time.dt,time.Tvec(k),model);
    
    [mX,PX]=MeanCov(Xquad,wquad);
    disp(['cond = ',num2str(cond(PX))])
    fullpdf=get_interpolated_pdf(X,probs,mX,PX,4);
    
    Xtestmc = mvnrnd(mX(:)',1^2*PX,1000);
    pXtest = fullpdf.func(Xtestmc);
    pX = fullpdf.func(X);
    logprobtestmc =log(pXtest);
    logprobX =log(pX);
    figure(10)
    plot3(X(:,1),X(:,2),log(probs),'ro',X(:,1),X(:,2),logprobX,'b+',Xtestmc(:,1),Xtestmc(:,2),logprobtestmc,'gs')
    hold on
    plot_nsigellip(mX(1:2),PX(1:2,1:2),1,'r',2)
    title(['k = ',num2str(k)])
    hold off
%     saveas(gcf,['nonorm_by4_gh5_k_',num2str(k)],'png')
    figure(11)
    plot3(X(:,1),X(:,2),probs,'ro',X(:,1),X(:,2),pX,'b+',Xtestmc(:,1),Xtestmc(:,2),pXtest,'gs')
    title(['k = ',num2str(k)])

    
    
%     [mX,PX]=MeanCov(Xquad,wquad);
%     [Xx,Xy]=meshgrid(linspace(mX(1)-2*sqrt(PX(1,1)),mX(1)+2*sqrt(PX(1,1)),25),linspace(mX(2)-2*sqrt(PX(2,2)),mX(2)+2*sqrt(PX(2,2)),25) );
%     Xp=[reshape(Xx,625,1),reshape(Xy,625,1)];
%     margprobs = fullpdf.func([Xp,repmat(mX(3:4)',625,1)]);
%     margprobs = get_2Dmarginalized_probs(Xp,1,2,X,probs,mX,PX,fullpdf,'dummyMC');
%     margprobs=reshape(margprobs,25,25);
%     
%     figure(1)
%     contour(Xx,Xy,margprobs,15)
%     title(['k = ',num2str(k)])
%     hold on
%     plot(XMC(:,1,k),XMC(:,2,k),'ro')
%     plot(X(:,1),X(:,2),'b*')
%     axis equal
%     axis square
%     hold off
%     
%     figure(2)
%     surf(Xx,Xy,margprobs)
%     title(['k = ',num2str(k)])
%     hold on
%     plot(XMC(:,1,k),XMC(:,2,k),'ro')
%     plot(X(:,1),X(:,2),'b*')
%     alpha 0.5
%     axis equal
%     axis square
%     hold off
    
    
    histXprior{k,1}=X;
    histXprior{k,2}=probs;
    
    pause(1)
    
    zk = model.h(Xtruth(k,:)')+sqrtm(model.R)*randn(model.hn,1);
    zk
    
    
    % do measurement update
    if k>2
        if rem(k,meas_freq_steps)==0
            disp("doing meas update")
            
            % %%%%%%%%%%%%%%%%% MEAS UPDATE %%%%%%%%%%%%%%%%%%%%%%
            [X,probs]=MeasUpdt_character(X,probs,zk,model);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(k),model);
            fullpdf=get_interpolated_pdf(X,probs,4);
            
            [mX,PX]=MeanCov(X(:,1:2),probs/sum(probs));
            [Xx,Xy]=meshgrid(linspace(mX(1)-9*sqrt(PX(1,1)),mX(1)+9*sqrt(PX(1,1)),50),linspace(mX(2)-9*sqrt(PX(2,2)),mX(2)+9*sqrt(PX(2,2)),50) );
            Xp=[reshape(Xx,625,1),reshape(Xy,625,1)];
            margprobs = get_2Dmarginalized_probs(Xp,1,2,X,probs,fullpdf,'dummyMC');
            margprobs=reshape(margprobs,25,25);
            
            figure(3)
            contour(Xx,Xy,margprobs)
            title(['k = ',num2str(k)])
            hold on
            plot(X(:,1),X(:,2),'b*')
            axis equal
            axis square
            hold off
            
            figure(4)
            surf(Xx,Xy,margprobs)
            title(['k = ',num2str(k)])
            hold on
            plot(X(:,1),X(:,2),'b*')
            alpha 0.5
            axis equal
            axis square
            hold off
            
        end
    end
    
    histXpost{k,1}=X;
    histXpost{k,2}=probs;
%     if k==12

keyboard
%     end
    
end