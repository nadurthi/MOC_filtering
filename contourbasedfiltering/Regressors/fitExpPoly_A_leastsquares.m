function mxentpoly = fitExpPoly_A_leastsquares(X,p,Pf,Xineq,XtestoutsideHull)
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
    Nextra = 100;
    extraFac = 0.1;
    size_top = 500;
end
if dim==6
    Nextra = 100;
    extraFac = 0.1;
    size_top = 1000;
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

Beqtrain = beq;

Aeqtrain=zeros(size(X,1),lamdim);
Aextra=zeros(size(Xextra,1),lamdim);
for ib=1:length(Pf)
    Aeqtrain(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},X);
    Aextra(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xextra);
end

Ntop = 100;
AeqTop = Aeqtrain(1:Ntop,:);
beqTop = Beqtrain(1:Ntop);

Xtest=[XtestoutsideHull;X;Xineq];
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
    Nineq=1;
    Aineq=zeros(1,lamdim);
    bineq=zeros(1,lamdim);
else
    Nineq = size(Xineq,1);
    Aineq=zeros(size(Xineq,1),lamdim);
    tic
    for ib=1:length(Pf)
        Aineq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xineq);
    end
    toc
    cc=min(beq);
    if cc<=0
        bineq = 1*cc*ones(size(Aineq,1),1);
    else
        bineq = 1*(cc+0)*ones(size(Aineq,1),1);
    end
    
end




%%
CC=linspace(0,10000,100);
Ntops=1:2:10;

%     CC=[0.5];
LAMS=zeros(lamdim,length(CC));
costs = zeros(1,length(CC));

% keyboard
Neqtrain = size(Aeqtrain,1);
brkflag=0;
cnt=0;

AA=[Aeqtrain;Aextra];  %Aineq
BB=[Beqtrain;bextra];  %bineq

ineqbnd=0.7*min(beq)+0.3*max(beq);
Bineqbnd = ineqbnd*ones(size(Aineq,1),1);
for ci = 1:length(CC)


    cnt = cnt + 1;
    
    lam = (AA'*AA+CC(ci)*eye(lamdim))\(AA'*BB);
  
    LAMS(:,ci)=lam;
    costs(ci) = norm(Aeqtrain*lam-Beqtrain,2);
    %         break
    
    [sum(Atest*lam <= 1.05*btest),sum(Aineq*lam<=Bineqbnd)]
    
    if isempty(Atest)==0
        if all(Atest*lam <= 1.05*btest) && all(Aineq*lam<=Bineqbnd) 
            disp('All probs are in the bounds')
            brkflag=1;
        end
    end
    figure(101)
    plot3(X(:,1),X(:,2),exp(beq),'ro')
    hold on
    plot3(X(:,1),X(:,2),exp(Aeqtrain*lam),'b+')
    hold off

    pause(1)

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

disp('done least poly fit')

end
