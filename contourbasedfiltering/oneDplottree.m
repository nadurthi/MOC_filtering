clc
clear
close all

mu = [1;5];
sigma = [1 1]; % shared diagonal covariance matrix

% gm = gmdistribution(mu,sigma);
f=@(x)0.3*normpdf(x,mu(1),sigma(1)^2)+0.7*normpdf(x,mu(2),sigma(2)^2);

a=-5;
b=12;

xplot = linspace(a,b,100);
yplot=f(xplot);

figure(1)
plot(xplot,yplot)

x=[];
cnt=1;
ErrRec=[];
while 1
    
    if isempty(x)
        xnew=a+(b-a)*rand(20,1);
    else
        xnew=[];
        nids=pmf_sampler(RMSEbox,1*5);
        for i=1:nb
                ns = sum(nids==i);
                lb = boxes{i,1};
                ub = boxes{i,2};
                xnew=[xnew;lb+(ub-lb)*rand(2*ns,1)];
        end
    end
    length(xnew)
    
    x=[x;xnew];
    y=f(x);
    
    tree = fitrtree(x,y,'MinParentSize',5,'MaxNumSplits',70,'MinLeafSize',5);
    ypred=tree.predict(x);
    boxes=getTree2Boxes(tree,a,b);
    
    nb=size(boxes,1);
    Mbox=zeros(1,nb);
    RMSEbox=zeros(1,nb);
    RelRMSEbox=zeros(1,nb);
    Nptbox=zeros(1,nb);
    A=0;
    for j=1:size(boxes,1)
        xlb=boxes{j,1};
        xub=boxes{j,2};
        m=boxes{j,4};
        A=A+(xub-xlb)*m;
        Mbox(j)=m;
        indin = (x>=xlb) & (x<xub);
        Nptbox(j)=sum(indin);
        RMSEbox(j)=sqrt(mean((ypred(indin)-y(indin)).^2));
        RelRMSEbox(j)=sqrt(mean(((ypred(indin)-y(indin))./(y(indin)+1)).^2));
    end
    Nptbox;
    figure(2)
    clf
    plot(xplot,yplot,'b','linewidth',2)
    hold on
    plot(x,y,'bo')
    disp(['Area = ',num2str(A)])
    for j=1:size(boxes,1)
        xlb=boxes{j,1};
        xub=boxes{j,2};
        yf=boxes{j,4};
        plot([xlb,xub],[yf,yf]/A,'g','linewidth',2)
        fill([xlb,xlb,xub,xub],[0,yf,yf,0]/A,'g','FaceColor','g','EdgeAlpha','0.1','FaceAlpha',0.4,'LineWidth',0.2)
        alpha 0.5
    end
    hold off
    grid on
    xlabel('x')
    ylabel('p(x)')
    plot_paper_optim(struct('setfigsize',1))
    pause(1)
    fname = ['simulations/TreeAdaptIllustrate/Approx_',num2str(cnt,'%03d')];
    saveas(gcf,fname,'png')
    saveas(gcf,fname,'fig')
    pause(1)
    
    RMSEerr=sqrt(mean((ypred-y).^2));
    RelRMSEerr=sqrt(mean(((ypred-y)./(y+1)).^2));
    
    ee=[cnt,length(x),RMSEerr,RelRMSEerr*100];
    ee
    ErrRec=vertcat(ErrRec,ee);
    
    cnt=cnt+1;
    if cnt>20
        break
    end
%     keyboard

    
end

figure(4)
yyaxis left
plot(ErrRec(:,1),ErrRec(:,3),'bo-','linewidth',2)
xlabel('iterations')
ylabel('RMSE')

yyaxis right
plot(ErrRec(:,1),ErrRec(:,2),'ro-','linewidth',2);
ylabel('# samples')

plot_paper_optim(struct('setfigsize',1))
grid
