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

model.h=@(x)[x(1);x(2)];
model.hn=2;
% model.R=diag([(0.1/constants.Re)^2,(0.5*pi/180)^2]);
model.R=diag([1^2,1^2]);
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
xf0=mvnrnd(x0(:)',P0);
% xf0 = x0;
Pf0 = P0;

Npf = 5000; %paricle filter points


% generate points on contours for characterisitc solutions


[X,w] = GH_points(zeros(model.fn,1),0.5^2*eye(model.fn),11);

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
% close all

X = Xinitial;
probs = probsinitial;

Xquad=Xquad_initial;
wquad=wquad_initial;

meas_freq_steps = 1;

histXprior=cell(length(time.Ntsteps),5);
histXpost=cell(length(time.Ntsteps),5);

histXprior{1,1} = X;
histXprior{1,2} = probs;

histXpost{1,1} = X;
histXpost{1,2} = probs;
teststeps = [33];

for k=2:time.Ntsteps
    k
    Xmctest = zeros(size(XMC,1),model.fn);
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
    disp([' k = ',num2str(k)])
    
    [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(k),model);
    [Xquad,wquad]=propagate_character(Xquad,wquad,time.dt,time.Tvec(k),model);
        
    [mX,PX]=MeanCov(X,probs/sum(probs));
    
%     [mX,PX]=MeanCov(Xquad,wquad);
%     mX=X(probs==max(probs),:);
%     PX=0;
%     wts=probs/sum(probs);
%     for i=1:size(X,1)
%        PX=PX+ wts(i)*(X(i,:)-mX)'*(X(i,:)-mX);
%     end
    disp(['cond = ',num2str(cond(PX))])
    fullnormpdf=NaN;
    
%         if any(k==teststeps)
    plotfolder='duffsim1';
    mkdir(plotfolder)
    plotsconf.plotfolder=plotfolder;
    plotsconf.nametag='prior';
    plotsconf.fig3.holdon = false;
    plotsconf.fig4.holdon = false;
    plotsconf.fig3.plottruth = true;
    plotsconf.fig4.plottruth = true;
    plotsconf.fig1.plottruth = true;
    plotsconf.fig2.plottruth = true;
    plotsconf.fig3.plotmeas = [];
    plotsconf.fig4.plotmeas = [];
    plotsconf.fig4.surfcol = 'green';
    plotsconf.fig3.contourZshift = -0.2;
    
    fullnormpdf=get_interp_pdf_0I_2D(X,probs,mX,PX,4,k,[],Xtruth(k,:),plotsconf); %Xtruth(k,:)
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
            [X,probs,Xquad,wquad,fullnormpdf]=MeasUpdt_character_2D(fullnormpdf,X,probs,Xquad,wquad,4,k,zk,Xtruth,model,Xmctest);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            histXpost{k,1}=X;
            histXpost{k,2}=probs;
            histXpost{k,3}=fullnormpdf;

            histXpost{k,4}=Xquad;
            histXpost{k,5}=wquad;
    
            
        end
    end
    

    

%     figure(49)
%     y=fullnormpdf.trueX2normX(X);
%     py=fullnormpdf.func(y);
%     probsXest=fullnormpdf.normprob2trueprob(py);
%     
%     plot3(X(:,1),X(:,2),probs,'ro',X(:,1),X(:,2),probsXest,'b+')
%     
%         if any(k==teststeps)
    
    keyboard
%         end
    
    
end

% save('sim1.mat')