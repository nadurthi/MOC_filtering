% test sparse grid on extrema of chebychev points
clc
close all
clear
d=6;
dim=d;
q=12;
pts1D=cell(1,q-d+2);
interpPoly1D = cell(1,q-d+2);

digits(50)



for i=1:q-d+2
    if i==1
%         [X,~] = GH_points(0,1,i);
%          [X,~] = GLeg_pts(i, -1, 1);
%          [ X, w1, w2 ] = kronrod ( i, 1e-12 );
%          X=[-X(1:end-1),X(end:-1:1)];
%         pts1D{i}=X(:);
% [ x, w ] = clenshaw_curtis_rule ( n, a, b, filename )
        pts1D{i}=chebychevextrema(1);
    else
%         [X,~] = GH_points(0,0.5^2,i);
%         pts1D{i}=X(:);
%          [X,~] = GLeg_pts(i, -1, 1);
%             [ X, w1, w2 ] = kronrod ( i, 1e-12 );
%          X=[-X(1:end-1),X(end:-1:1)];
%          
%          pts1D{i}=X(:);
        
        pts1D{i}=chebychevextrema(2^(i-1)+1);
%         pts1D{i}=chebychevextrema(i);
    end
    tic
%     interpPoly1D{i}=LagrangePoly1D( pts1D{i} ,1,1);  % vpa(pts1D{i})  pts1D{i}
    interpPoly1D{i}=LagrangePoly1D_matpoly( pts1D{i} ,1,1);  % vpa(pts1D{i})  pts1D{i}
    
    toc
end
disp('done')

f=@(x)sin(x);
X=pts1D{end};
c=f(X);
P=LinearComb_VectorOfPolys(c,interpPoly1D{end});
[c,evaluate_polyND(P,X),c-evaluate_polyND(P,X)]

X=sparseGridConstructor(q,d,pts1D);
size(X)

plot(X(:,1),X(:,2),'bo')

%%
f=@(x)0.25*mvnpdf(x,[0.2,0.2],0.2^2*eye(2))+0.25*mvnpdf(x,-[0.5,0.5],0.2^2*eye(2))+0.25*mvnpdf(x,[-0.5,0.5],0.2^2*eye(2))+0.25*mvnpdf(x,[0.5,-0.5],0.2^2*eye(2));
% f=@(x)ones(size(x,1),1);
FuncTable = [X,f(X)];
FuncTablelog = [X,log(f(X))];

figure
plot3(X(:,1),X(:,2),FuncTable(:,3),'bo')

PnD=sparseProductInterpPoly(q,d,pts1D,interpPoly1D,FuncTablelog);

Festlog=evaluate_polyND(PnD,X);
[Festlog,FuncTablelog(:,3)]
PnD



[xx,yy]=meshgrid(linspace(-1,1,105));
Fplotlog = zeros(size(xx));
Ftrue = zeros(size(xx));
for i=1:size(xx,1)
    for j=1:size(xx,2)
        Fplotlog(i,j) = evaluate_polyND(PnD,[xx(i,j),yy(i,j)]);
        Ftrue(i,j) = f([xx(i,j),yy(i,j)]);
    end
end

Fplot = exp(Fplotlog);
% Fplot = (Fplotlog);

figure
plot3(X(:,1),X(:,2),FuncTable(:,3),'bo')
hold on
surf(xx,yy,Fplot,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
camlight right; lighting phong
alpha 0.5

surf(xx,yy,Ftrue,'FaceColor','g','EdgeColor','none','FaceAlpha',0.5);
camlight right; lighting phong
alpha 0.5


%%
close all
Pf=Basis_polyND(dim,15);
lamdim = length(Pf);
Aeq=zeros(size(X,1),lamdim);
for ib=1:length(Pf)
    Aeq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},X);
end

beq=FuncTablelog(:,3);
Neq=length(beq);
wts=FuncTable(:,3)/sum(FuncTable(:,3));

cvx_begin
variables lam(lamdim)  %teqTop(Ntop)
variable t
minimize( 1*norm(lam,1) +1000*abs(t) ) %+500*norm(teqTop,2)+500*norm(teq,2)
subject to
Aeq*lam-beq<=t;
Aeq*lam-beq>=-t;
cvx_end

mxentpoly=zeros(1,dim+1);
for i=1:lamdim
    mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lam(i),Pf{i}),'add');
end
mxentpoly=simplify_polyND(mxentpoly);

[evaluate_polyND(mxentpoly,X),FuncTablelog(:,3)]


Fcvxlog = zeros(size(xx));
for i=1:size(xx,1)
    for j=1:size(xx,2)
        Fcvxlog(i,j) = evaluate_polyND(mxentpoly,[xx(i,j),yy(i,j)]);
    end
end

Fcvx = exp(Fcvxlog);
% Fplot = (Fplotlog);

figure
plot3(X(:,1),X(:,2),FuncTable(:,3),'bo')
hold on
surf(xx,yy,Fplot,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
camlight right; lighting phong
alpha 0.5

surf(xx,yy,Ftrue,'FaceColor','g','EdgeColor','none','FaceAlpha',0.5);
camlight right; lighting phong
alpha 0.5

surf(xx,yy,Fcvx,'FaceColor','b','EdgeColor','none','FaceAlpha',0.5);
camlight right; lighting phong
alpha 0.5
