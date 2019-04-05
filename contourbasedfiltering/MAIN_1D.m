%%
% Satellite problem

clc
close all
clear

format longg

digits(50)

%% time


time.t0=0;
time.tf=20;
time.dt=1;

time.Tvec=time.t0:time.dt:time.tf;
time.Ntsteps=length(time.Tvec);

%% models

model.f=@(dt,tk,xk)xk/2+25*xk./(1+xk.^2)+8*cos(1.2*tk);
model.fn=1;
model.Q = diag([(1e-6)^2]);
model.fstates = {'x'};

% model.h=@(x)[sqrt(x(1)^2+x(2)^2)];
model.h=@(x)x.^2/20;
model.hn=1;
model.hvec=@(x)model.h(x(:));

% model.hstates = {'x','y'};
% model.hstates = {'r','th'};
model.hstates = {'z'};
% model.R=diag([(0.1/constants.Re)^2,(0.5*pi/180)^2]);
% model.R=diag([0.1^2,0.1^2]);
% model.R=diag([1^2,(2*pi/180)^2]);
model.R=diag([(1)^2]);
% model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);


%% generate truth

x0=0;
P0=5^2*eye(1);

Xtruth = zeros(time.Ntsteps,model.fn);
Xtruth(1,:)=x0;
for k=2:time.Ntsteps
    Xtruth(k,:)=model.f(time.dt,time.Tvec(k),Xtruth(k-1,:));
end

Ntsteps=1000;
plotTvec=linspace(time.t0,time.tf,Ntsteps);
plotdt = plotTvec(2)-plotTvec(1);
Xtruplot = zeros(Ntsteps,model.fn);
Xtruplot(1,:)=x0;
for k=2:Ntsteps
    Xtruplot(k,:)=model.f(plotdt,plotTvec(k),Xtruplot(k-1,:));
end


figure
plot(time.Tvec,Xtruth(:,1),'ro-')
axis equal
axis square



% plotting the propagatin of MC
Nmc=10000;
XMC=zeros(Nmc,model.fn,time.Ntsteps);
XMC(:,:,1)=mvnrnd(x0',P0,Nmc);
for i=1:Nmc
    for k=2:time.Ntsteps
        XMC(i,:,k)=model.f(time.dt,time.Tvec(k),XMC(i,:,k-1));
    end
end

figure
for ti=1:time.Ntsteps
    hist(XMC(:,1,ti),linspace(-25,25,200))
    xlim([-25,25])
    pause(1)
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

GMMinitial=getInitialGMM(xf0(:),0.7^2*Pf0,7);
model.gmmmethod='ut';

% generate points on contours for characterisitc solutions

model.pointGenerator = @(mx,Px)GH_points(mx,Px,11);


% model.pointGenerator = @(mx,Px)mvnrnd(mx,Px,200);
model.pointscalefactor = 0.4;
[X,w] = model.pointGenerator(zeros(model.fn,1),model.pointscalefactor^2*eye(model.fn));
% [X,w] = mvnrnd(zeros(model.fn,1),1^2*eye(model.fn),100);

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

GMM = GMMinitial;

meas_freq_steps = 1;

histXprior=cell(time.Ntsteps,5);
histXpost=cell(time.Ntsteps,5);

histXprior{1,1} = X;
histXprior{1,2} = probs;

histXpost{1,1} = X;
histXpost{1,2} = probs;
teststeps = [33];

plotfolder='simulations/duffsim_thmeas_test';
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
%     [Xquad,wquad]=propagate_character(Xquad,wquad,time.dt,time.Tvec(k),model);
    [xfquad,Pfquad]=QuadProp(xfquad,Pfquad,time.dt,time.Tvec(k),model,'ut');

    GMM=prior_prop_GMM(GMM,time.dt,time.Tvec(k),model);
    

    
    

    disp(['cond = ',num2str(cond(PX))])

    [mX,PX]=MeanCov(X,probs/sum(probs));
    fullnormpdf=get_interp_pdf_0I_duff(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)
%     fullnormpdf=get_interp_pdf_0I_2D(X,probs,mX,PX,4,k,[],Xtruth(k,:),plotsconf); %Xtruth(k,:)
    priorpdfnorm=fullnormpdf;
    
    plotpdfs_prior_2D(k,fullnormpdf,X,probs,xfquad,Pfquad,GMM,Xmctest,Xtruth(k,:),savePriorProps)
    
    %     keyboard
    
%     figure(12)
%     hold on
%     [mG,PG]=GMMcell2array(GMM);
%     plot(mG(:,1),mG(:,2),'b*','MarkerSize',10)
%     hold off
    
%     X=histXprior{k,1};
%     probs=histXprior{k,2};
%     fullnormpdf=histXprior{k,3};
    
    
    histXprior{k,1}=X;
    histXprior{k,2}=probs;
    histXprior{k,3}=fullnormpdf;
    
    histXprior{k,4}=xfquad;
    histXprior{k,5}=Pfquad;
    
   
    pause(1)
    
    
%     X=histXprior{k,1};
%     probs=histXprior{k,2};
    
    % do measurement update
    if k>=2
        if rem(k,meas_freq_steps)==0
            disp("doing meas update")
            zk = model.h(Xtruth(k,:)')+sqrtm(model.R)*randn(model.hn,1);
            zk
            
            % %%%%%%%%%%%%%%%%% MEAS UPDATE %%%%%%%%%%%%%%%%%%%%%%
            [X,probs,fullnormpdf]=MeasUpdt_character_modf_duff(X,probs,priorpdfnorm,4,k,zk,Xtruth(k,:),model,Xmctest,11);
%             [X,probs,fullnormpdf]=MeasUpdt_character_modf(fullnormpdf,X,probs,4,k,zk,Xtruth,model,Xmctest);
            [xfquad,Pfquad]=QuadMeasUpdt(xfquad,Pfquad,zk,time.dt,time.Tvec(k),model,'ut');
            
            GMM=post_GMM(GMM,zk,time.dt,time.Tvec(k),model);
            
            plotpdfs_post_2D(k,priorpdfnorm,fullnormpdf,model,Xtruplot,X,probs,xfquad,Pfquad,GMM,Xmctest,Xtruth(k,:),zk,savePostProps)
%             plotpdfs_post_2D(Tk,pdfnormprior,pdfnorm,model,XtruthAll,X,probs,mquadf,Pquadf,GMM,Xmc,Xtruth,zk,saveprops)
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            histXpost{k,1}=X;
            histXpost{k,2}=probs;
            histXpost{k,3}=fullnormpdf;

            histXpost{k,4}=xfquad;
            histXpost{k,5}=Pfquad;
            histXpost{k,6}=zk;
            
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