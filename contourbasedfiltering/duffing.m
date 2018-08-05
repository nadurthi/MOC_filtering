clc
clear all
close all

mu0=[5,5];
P0=0.7^2*eye(2);
[t,ground]=ode45(@duff,linspace(0,10,1000),mu0);


model.Ngrid=101;

[x1,x2]=meshgrid(linspace(mu0(1)-2,mu0(1)+2,model.Ngrid),linspace(mu0(2)-6,mu0(2)+6,model.Ngrid));
ps=zeros(size(x1));
[a,b]=size(x1);
model.a=a;
model.b=b;

for i=1:size(x1,1)
    for j=1:size(x1,2)
        ps(i,j) = mvnpdf([x1(i,j),x2(i,j)],mu0,P0);
    end
end

X=[reshape(x1, a*b,1 ),reshape(x2, a*b,1 )];
probs = reshape(ps, a*b,1 );

model.f = @duff_prop_model;
model.R=diag([0.5^2,1^2]);
model.h=@(x)[x(1);x(2)];
model.hn=2;
model.z_pdf = @(z,x)mvnpdf(z,model.h(x),model.R);

dt=0.5;

model.dt=dt;
model.t0=0;
model.tf=10;
model.timesteps=model.t0:model.dt:model.tf;

[t,Xtruth]=ode45(@duff,model.timesteps,mu0);

% figure
% plot(ground(1:200,1),ground(1:200,2))

%%

% for k=1:10
%     [X,probs]=propagate_character(X,probs,dt,model);
%     
%     figure
%     hold on
%     mesh(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b) )
%     plot(ground(:,1),ground(:,2))
% end
%%
% running the filter

meas_freq_steps = 3;
histXprior=cell(length(model.timesteps),2);
histXpost=cell(length(model.timesteps),2);

histXprior{1,1} = X;
histXprior{1,2} = probs;

histXpost{1,1} = X;
histXpost{1,2} = probs;


for k=2:length(model.timesteps)
    k
    
    
    [X,probs]=propagate_character(X,probs,dt,model);
    
    histXprior{k,1}=X;
    histXprior{k,2}=probs;
    
    figure(1)
    hold off
    % axis equal
    mesh(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b) )   
    hold on
    for i=1:k
        mesh(reshape(histXprior{i,1}(:,1),a,b),reshape(histXprior{i,1}(:,2),a,b),reshape(histXprior{i,2},a,b) )
    end

    plot(ground(1:200,1),ground(1:200,2))
    title('prior')
    xlabel('x_1')
    ylabel('x_2')
    plot_prop_paper
    
    figure(21)
    hold off
    % axis equal
    surf(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b),'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
    camlight right; lighting phong  
    
    hold on
    for i=1:k
        surf(reshape(histXprior{i,1}(:,1),a,b),reshape(histXprior{i,1}(:,2),a,b),reshape(histXprior{i,2},a,b),'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
        camlight right; lighting phong
    end

    plot(ground(1:200,1),ground(1:200,2))
    title('prior')
    xlabel('x_1')
    ylabel('x_2')
    plot_prop_paper
    

    
    % axis equal
    
%     figure(1)
%     hold off
%     % axis equal
%     mesh(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b) )   
%     hold on
%     for i=1:k
%         mesh(reshape(histXprior{i,1}(:,1),a,b),reshape(histXprior{i,1}(:,2),a,b),reshape(histXprior{i,2},a,b) )
%     end
% 
%     plot(ground(1:200,1),ground(1:200,2))
%     title('prior')
%     xlabel('x_1')
%     ylabel('x_2')
%     plot_prop_paper
    
    
%     plot(X(:,1),X(:,2),'ko')


    pause(1)
    % generate measurement
    zk = model.h(Xtruth(k,:)')+sqrtm(model.R)*randn(model.hn,1);
    zk
    
    
    
    % contour plots
    figure(3)
    hold off
%     axis equal
    contour(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b),10 )
    hold on
    for i=1:k
        contour(reshape(histXprior{i,1}(:,1),a,b),reshape(histXprior{i,1}(:,2),a,b),reshape(histXprior{i,2},a,b),10 )
    end
    plot(ground(1:200,1),ground(1:200,2))
    plot_prop_paper
    title('prior contour')
    xlabel('x_1')
    ylabel('x_2')
    plot_prop_paper
%     axis equal
    
    % do measurement update
    if k>2
        if rem(k,meas_freq_steps)==0
            disp("doing meas update")
            
            figure(9)
            hold off
            surf(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b),'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
            camlight right; lighting phong
            hold on
            title('prior and posterior')
            xlabel('x_1')
            ylabel('x_2')
            plot_prop_paper
            
            figure(29)
            hold off
            mesh(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b));
            hold on
            title('prior and posterior')
            xlabel('x_1')
            ylabel('x_2')
            plot_prop_paper
            
            figure(39)
            hold off
            mesh(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b));
            colormap summer
            hold on
            for i=1:k-1
                mesh(reshape(histXprior{i,1}(:,1),a,b),reshape(histXprior{i,1}(:,2),a,b),reshape(histXprior{i,2},a,b) )
                colormap summer
            end

            plot(ground(1:200,1),ground(1:200,2))

            title('prior and posterior')
            xlabel('x_1')
            ylabel('x_2')
            plot_prop_paper
            
            figure(49)
            hold off
            surf(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b),'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
            camlight right; lighting phong
            hold on
            for i=1:k-1
                surf(reshape(histXprior{i,1}(:,1),a,b),reshape(histXprior{i,1}(:,2),a,b),reshape(histXprior{i,2},a,b),'FaceColor','green','EdgeColor','none','FaceAlpha',0.7);
                camlight right; lighting phong
            end

            plot(ground(1:200,1),ground(1:200,2))

            title('prior and posterior')
            xlabel('x_1')
            ylabel('x_2')
            plot_prop_paper
            
            
            figure(4)
            hold off
            mesh(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b) )
            hold on
            plot(zk(1),zk(2),'k*','MarkerSize',6)
            [Z1,Z2]=meshgrid(linspace(zk(1)-2,zk(1)+2,101),linspace(zk(2)-5,zk(2)+5,101));
            pZ=zeros(size(Z1));
            for i=1:size(Z1,1)
                for j=1:size(Z1,2)
                    pZ(i,j)=mvnpdf([Z1(i,j),Z2(i,j)], zk(:)', model.R);
                end
            end
%             [Zin1,Zin2]=meshgrid(linspace(Xtruth(k,1)-2,Xtruth(k,1)+2,101),linspace(Xtruth(k,2)-5,Xtruth(k,2)+5,101));
            surf(Z1,Z2,pZ,'FaceAlpha',0.8,'EdgeColor','flat')
            title('prior and likelihood')
            xlabel('x_1')
            ylabel('x_2')
            plot_prop_paper
            
            figure(7)
            hold off
            contour(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b),10 )
            hold on
            plot(zk(1),zk(2),'k*','MarkerSize',6)
            contour(Z1,Z2,pZ,10)
            title('prior and likelihood contours')
            xlabel('x_1')
            ylabel('x_2')
            plot_prop_paper
            
            % %%%%%%%%%%%%%%%%% MEAS UPDATE %%%%%%%%%%%%%%%%%%%%%%
            [X,probs]=MeasUpdt_character(X,probs,zk,model);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            figure(9)
            surf(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b),'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
            camlight right; lighting phong
            
            figure(29)
            mesh(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b))
            
            figure(39)
            mesh(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b))
            colormap winter
            
            figure(49)
            surf(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b),'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
            camlight right; lighting phong
            
            
            
            figure(2)
            hold off
            mesh(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b) )
            hold on
            plot(zk(1),zk(2),'k*','MarkerSize',6)
            plot(ground(1:200,1),ground(1:200,2))
            title('posterior')
            xlabel('x_1')
            ylabel('x_2')
            plot_prop_paper

            figure(5)
            hold off
            contour(reshape(X(:,1),a,b),reshape(X(:,2),a,b),reshape(probs,a,b),10 )
            hold on
            plot(zk(1),zk(2),'k*','MarkerSize',6)
            plot(ground(1:200,1),ground(1:200,2))
            title('posterior contours')
            xlabel('x_1')
            ylabel('x_2')
            plot_prop_paper
            
            
            
        end
    end
    
    histXpost{k,1}=X;
    histXpost{k,2}=probs;
    
    keyboard
    
end

