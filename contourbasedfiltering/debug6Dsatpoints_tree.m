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
time.dt=0.5*60*60 *constants.trueT2normT;
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
model.fstates = {'a','e','i','Om','om','M'};
% model.fstates = {'x','y','z','vx','vy','vz'};

model.h=@(x)radmodel(xk,1);
model.hn=2;
model.hvec=@(x)sensmodel2vec(x,model.h,model.hn);
model.R=diag([(0.1*pi/180)^2,(0.1*pi/180)^2]);
model.hstates = {'azi','el'};
% model.R=diag([(0.1/constants.Re)^2]);
model.z_pdf =  @(z,x)mvnpdf(z,model.h(x),model.R);



%% generate truth

x0=[10000,10,3000,-1,7.5,3.5]';
OE = cart2orbelem(x0,constants.mu)
% [ r, v, Ehat ] = FnG(0, time.dt, x0(1:3), x0(4:6), 1);
% [ r1, v1, Ehat ] = FnG(0, time.dt, r, -v, 1);

P0=diag([0.001^2,0.001^2,0.001^2,0.00001^2,0.00001^2,0.00001^2]);

x0(1:3)=x0(1:3)*constants.trueX2normX;
x0(4:6)=x0(4:6)*constants.trueV2normV;
P0(1:3,1:3)=P0(1:3,1:3)*constants.trueX2normX^2;
P0(4:6,4:6)=P0(4:6,4:6)*constants.trueV2normV^2;

X=mvnrnd(x0(:)',P0,1000);
nn=size(X,1);


Xtruth = zeros(time.Ntsteps,model.fn);
Xtruth(1,:)=x0;
for k=2:time.Ntsteps
    Xtruth(k,:)=model.f(time.dt,time.Tvec(k-1),Xtruth(k-1,:));
end




% plotting the propagatin of MC
Nmc=100000;
XMC=zeros(Nmc,model.fn,time.Ntsteps);
% XMCcoe=zeros(Nmc,model.fn,time.Ntsteps);
XMC(:,:,1)=mvnrnd(x0',P0,Nmc);
% for i=1:Nmc
% %     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe(XMC(i,1:3,1),XMC(i,4:6,1), 1);
% %     XMCcoe(i,:,1)=[p,a,ecc,incl,omega,argp];
% end
for i=1:Nmc
    i
    for k=2:time.Ntsteps
        XMC(i,:,k)=model.f(time.dt,time.Tvec(k-1),XMC(i,:,k-1));
    end
end
pMC = mvnpdf(XMC(:,:,1),x0',P0);



%% 
Xmctest=zeros(size(XMC,1),6);
for k=1:time.Ntsteps
    for ii=1:size(XMC,1)
        Xmctest(ii,:) = XMC(ii,:,k);
    end
    [m,p]=MeanCov(Xmctest,ones(Nmc,1)/Nmc);
    A=sqrtm(inv(p));

    for ii=1:size(XMC,1)
        Xmctest(ii,:) = A*(Xmctest(ii,:)-m(:)')';
    end
    
    figure(1)
    plot(Xmctest(:,1)*constants.normX2trueX,Xmctest(:,2)*constants.normX2trueX,'r.')
    pause(1)
    k
end
%%
figure
subplot(3,2,1)
plot(Xmctest(:,1),Xmctest(:,2),'r.')
xlabel('x')
ylabel('y')

subplot(3,2,3)
plot(Xmctest(:,2),Xmctest(:,3),'r.')
xlabel('y')
ylabel('z')


subplot(3,2,5)
plot(Xmctest(:,3),Xmctest(:,4),'r.')
xlabel('z')
ylabel('vx')

subplot(3,2,2)
plot(Xmctest(:,4),Xmctest(:,5),'r.')
xlabel('vx')
ylabel('vy')

subplot(3,2,4)
plot(Xmctest(:,5),Xmctest(:,6),'r.')
xlabel('vy')
ylabel('vz')

subplot(3,2,6)
plot(Xmctest(:,1),Xmctest(:,4),'r.')
xlabel('x')
ylabel('vx')

