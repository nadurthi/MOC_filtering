mx=[2,2];
Px=[3,1;1,3];
X=mvnrnd(mx,Px,100);
probs = mvnpdf(X,mx,Px);

figure
plot3(X(:,1),X(:,2),probs,'ro')
grid
axis([-3,6,-3,6,0,0.1])

%%
[m,Pcov] = MeanCov(X,probs/sum(probs));


%%
[pdf,pdf0Inorm] = get_interp_pdf(X,probs,4);


[probs,pdf.func(X)]

%% test normalization
copypdf =pdf;
copypdf.poly = update_or_insert_coeff_NDpoly(copypdf.poly,zeros(1,size(X,2)),5);
copypdf.func=@(x)exp(evaluate_polyND(copypdf.poly,x));

normpdf = normalize_exp_pdf(copypdf,X,'dummyMC')
normpdf.poly
% [probs,pdf.func(X),normpdf.func(X)]


%% test better interpolation and normalization
[mX,PX]=MeanCov(X,probs/sum(probs));
Y=zeros(size(X));
A=sqrtm(inv(PX));
for i=1:size(X,1)
   Y(i,:) = A*(X(i,:)'-mX) ;
end
figure
plot3(X(:,1),X(:,2),probs,'r+')


figure
plot3(Y(:,1),Y(:,2),probs/det(A),'bo')


%% 

ss=[0.775538187751073,1.32,-0.416809885467009,-0.183497676615338];
Xnn=mvnrnd(ss,0.01*PX);

Idx = knnsearch(Xtestmc,ss,'K',20);
Y=Xtestmc(Idx,:);


p1=3;p2=4;
figure
plot(Y(:,p1),Y(:,p2),'ro',ss(p1),ss(p2),'b^','MarkerSize',4)
