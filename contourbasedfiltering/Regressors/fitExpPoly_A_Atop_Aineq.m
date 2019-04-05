function mxentpoly = fitExpPoly_A_Atop_Aineq(X,p,Pf,Xineq,XtestoutsideHull)
[N,dim] = size(X);
lamdim = length(Pf);

% keyboard

factconst = max(p)/10;
pnfit = p/factconst;
beq = log(pnfit);

Nextra = dim;
if dim==2
    Nextra = 100;
    extraFac = 0.05;
    size_top = N;
end
if dim==4
    Nextra = 10;
    extraFac = 0.05;
    size_top = 1000;
end
if dim==6
    Nextra = 20;
    extraFac = 0.1;
    size_top = 200;
end
Xextra=zeros(Nextra*size_top,size(X,2));
bextra = zeros(Nextra*size_top,1);
for i=1:size_top
    Xextra((i-1)*Nextra+1:i*Nextra,:) = mvnrnd(X(i,:),extraFac^2*eye(dim),Nextra);
    %         bextra(i) = beq(i);
end
idx = knnsearch(X,Xextra,'K',3);
for j=1:size(Xextra,1)
    dists = sqrt(sum((X(idx(j,:),:)-repmat(Xextra(j,:),3,1)).^2,2));
    wts = 1./dists;
    wts = wts/sum(wts);
    bextra(j) = sum(wts.*beq(idx(j,:)));
end

%     Xtrain = [X;Xextra];
%     Xtrain=X;
%     Beqtrain = [beq;bextra];

Pgrads=cell(length(Pf),dim);
for ib=1:length(Pf)
    for id=1:dim
        Pgrads{id,ib}=diff_polyND(Pf{ib},id);
    end
end
% Pgrads=reshape(Pgrads,length(Pf)*dim,1);

Beqtrain = beq;

Aeqtrain=zeros(size(X,1),lamdim);
Aextra=zeros(size(Xextra,1),lamdim);

for ib=1:length(Pf)
    Aeqtrain(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},X);
    Aextra(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xextra);
end

Ntop = 40;
AeqTop = Aeqtrain(1:Ntop,:);
beqTop = Beqtrain(1:Ntop);

Xtest=[XtestoutsideHull;X];
if isempty(Xtest)
    Atest=[];
    btest=[];
else
    Atest=zeros(size(Xtest,1),lamdim);
    for ib=1:length(Pf)
        Atest(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xtest);
    end
    btest = max(beq)*ones(size(Atest,1),1);
end

if isempty(Xineq)
    Nineq = 1;
    Aineq = zeros(1,lamdim);
    bineq = zeros(1,lamdim);
    Adervs = zeros(1,lamdim);

else
    Nineq = size(Xineq,1);
    Aineq=zeros(size(Xineq,1),lamdim);
    Adervs=zeros(size(Xineq,1)*dim,lamdim);

    for ib=1:length(Pf)
        for id=1:dim
            Adervs((id-1)*Nineq+1:(id)*Nineq,ib)=evaluate_PolyforLargetSetX(Pgrads{id,ib},Xineq);
        end
    end
    
    
    tic
    for ib=1:length(Pf)
        Aineq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xineq);
    end
    toc
    cc=min(beq);
    cc
    if cc<=0
        bineq = -(10)*ones(size(Aineq,1),1);
    else
        bineq = (0.2*cc+0.8*(-10))*ones(size(Aineq,1),1);
    end
    
end




%%
CC=[linspace(0.1,1000,10),linspace(2000,15000,10)];
% Ntops=1:2:floor(lamdim/2);
Ntops=[1:10];

%     CC=[0.5];
LAMS=zeros(lamdim,length(CC));
costs = zeros(1,length(CC));

% keyboard
Neqtrain = size(Aeqtrain,1);
brkflag=0;
cnt=0;

dervsmax=10;

for Ntop=Ntops
    for ci = 1:length(CC)
    
%         if cnt> 200
%             keyboard
%         end
        cnt = cnt + 1;
        
        %         [lam,teq]=polyfitcvx_parasolve(CC(ci),lamdim,Aeqtrain,Beqtrain,AeqTop,beqTop,Aineq,bineq);
        [lam,teq]=polyfitcvx_parasolve_extra(CC(ci),lamdim,Aeqtrain,Beqtrain,AeqTop(1:Ntop,:),beqTop(1:Ntop),Aineq,bineq,Aextra,bextra);
%         [lam,teq]=polyfitcvx_parasolve_extra_dervconstraints(CC(ci),lamdim,Aeqtrain,Beqtrain,AeqTop(1:Ntop,:),beqTop(1:Ntop),Aineq,bineq,Aextra,bextra,Adervs,dervsmax)
        %     [lam,teq]=polyfitcvx_parasolve_extra(CC(ci),lamdim,Aeqtrain,Beqtrain,AeqTop,beqTop,[],[],Aextra(1:10,:),bextra(1:10));
        
        LAMS(:,ci)=lam;
        costs(ci) = norm(teq,2);
        %         break
        if isempty(Atest)==0
            if all(exp(Atest*lam) <= 1.5*exp(btest) )
                disp('All probs are in the bounds')
                brkflag=1;
            end
        end
        figure(201)
        plot3(X(:,1),X(:,2),exp(beq),'ro')
        hold on
        plot3(X(:,1),X(:,2),exp(Aeqtrain*lam),'b+')
        hold off
        title([num2str(Ntop),',',num2str(ci)])
        
        if dim>=4
            figure(202)
            plot3(X(:,3),X(:,4),exp(beq),'ro')
            hold on
            plot3(X(:,3),X(:,4),exp(Aeqtrain*lam),'b+')
            hold off
            title([num2str(Ntop),',',num2str(ci)])
        end
        
        pause(1)
        
        if brkflag==1
            break
        end
    end
    if brkflag==1
        break
    end
end

lam(1) = lam(1)+log(factconst);
lamsol = lam;

mxentpoly=zeros(1,dim+1);
for i=1:lamdim
    mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
end
mxentpoly=simplify_polyND(mxentpoly);

disp('done poly fit')

end
%% -------------------------------------------------------------------
% function [lam,teq]=polyfitcvx_parasolve(C,lamdim,Aeq,beq,AeqTop,beqTop,Aineq,bineq)
% Neq = size(Aeq,1);
% Ntop = size(AeqTop,1);
% 
% cvx_begin
% variables teq(Neq) lam(lamdim) teqTop(Ntop)
% minimize( C*norm(lam,1)+50*norm(teq,2)+100*norm(teqTop,2) )
% subject to
% Aeq*lam==beq+teq;
% AeqTop*lam==beqTop+teqTop;
% Aineq*lam <=bineq;
% cvx_end
% end

function [lam,teq]=polyfitcvx_parasolve_extra(C,lamdim,Aeq,beq,AeqTop,beqTop,Aineq,bineq,Aextra,bextra)
Nextra = size(Aextra,1);
Neq = size(Aeq,1);
Ntop = size(AeqTop,1);
Nineq = size(Aineq,1);
wts=exp(beq)+1;
if Nineq>0
    cvx_begin
    variables teq(Neq) lam(lamdim) textra(Nextra)  %  teqTop(Ntop)
    minimize( C*norm(lam,1)+10*norm(textra,2)+50*norm(teq./(wts),2))  % +500*norm(teqTop,2)
    subject to
    Aeq*lam==beq+teq;
    AeqTop*lam==beqTop; %+teqTop
    Aextra*lam == bextra+textra;
    Aineq*lam <=bineq;
    cvx_end
else
    cvx_begin
    variables teq(Neq) lam(lamdim)  textra(Nextra) %teqTop(Ntop)
    minimize( C*norm(lam,1)+10*norm(textra,2)+50*norm(teq,2) ) %+500*norm(teqTop,2)
    subject to
    Aeq*lam==beq+teq;
    AeqTop*lam==beqTop; %+teqTop
    Aextra*lam == bextra+textra;
    cvx_end
end
end
function [lam,teq]=polyfitcvx_parasolve_extra_dervconstraints(C,lamdim,Aeq,beq,AeqTop,beqTop,Aineq,bineq,Aextra,bextra,Adervs,dervsmax)
Nextra = size(Aextra,1);
Neq = size(Aeq,1);
Ntop = size(AeqTop,1);
Nineq = size(Aineq,1);
if Nineq>0
    cvx_begin
    variables teq(Neq) lam(lamdim) textra(Nextra) tderiv %  teqTop(Ntop)
    minimize( C*norm(lam,1)+10*norm(textra,2)+50*norm(teq,2)-tderiv)  % +500*norm(teqTop,2)
    subject to
    Aeq*lam==beq+teq;
    AeqTop*lam==beqTop; %+teqTop
    Aextra*lam == bextra+textra;
    Aineq*lam <=bineq;
    Adervs*lam <=-tderiv;
    tderiv>=-1;
    cvx_end
else
    cvx_begin
    variables teq(Neq) lam(lamdim)  textra(Nextra) %teqTop(Ntop)
    minimize( C*norm(lam,1)+10*norm(textra,2)+50*norm(teq,2) ) %+500*norm(teqTop,2)
    subject to
    Aeq*lam==beq+teq;
    AeqTop*lam==beqTop; %+teqTop
    Aextra*lam == bextra+textra;
    Adervs*lam <=dervsmax;
    cvx_end
end
end