[pdf,~] = get_interp_pdf(X,probs,4);
figure(10)
plot3(X(:,1),X(:,2),log(probs),'ro',X(:,1),X(:,2),log(pdf.func(X)),'b+')
figure(11)
plot3(X(:,1),X(:,2),probs,'ro',X(:,1),X(:,2),pdf.func(X),'b+')


fullpdf=get_interpolated_pdf(X,probs,4);
figure(10)
plot3(X(:,1),X(:,2),log(probs),'ro',X(:,1),X(:,2),log(fullpdf.func(X)),'b+')
figure(11)
plot3(X(:,1),X(:,2),probs,'ro',X(:,1),X(:,2),fullpdf.func(X),'b+')

[mX,PX]=MeanCov(X(:,1:2),probs/sum(probs));
[Xx,Xy]=meshgrid(linspace(mX(1)-9*sqrt(PX(1,1)),mX(1)+9*sqrt(PX(1,1)),25),linspace(mX(2)-9*sqrt(PX(2,2)),mX(2)+9*sqrt(PX(2,2)),25) );
Xp=[reshape(Xx,625,1),reshape(Xy,625,1)];
margprobs = get_2Dmarginalized_probs(Xp,1,2,X,probs,fullpdf,'dummyMC');
margprobs=reshape(margprobs,25,25);

figure(1)
contour(Xx,Xy,margprobs)
title(['k = ',num2str(k)])
hold on
plot(XMC(:,1,k),XMC(:,2,k),'ro')
plot(X(:,1),X(:,2),'b*')
axis equal
axis square
hold off

figure(2)
surf(Xx,Xy,margprobs)
title(['k = ',num2str(k)])
hold on
plot(XMC(:,1,k),XMC(:,2,k),'ro')
plot(X(:,1),X(:,2),'b*')
alpha 0.5
axis equal
axis square
hold off


%%

[pdf,pdf0Inorm] = get_interp_pdf_nonlin(X,probs,mquad,Pquad,4)