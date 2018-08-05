function f=integrate_func_exppdf(F,pdf,m,P,method)
% integrate F with respect to the exp pdf
% exp pdf has a rough mean and covariance as m and P
% m,P are used to normalize the exp pdf
% *** TODO:  mix of gauss integration for beter accurqacy
ndim = length(m);
m=m(:);

% given in p(x)
% y=A(x-mu)
% x = inv(A)*y+mu

A=sqrtm(inv(P));
Ainv = inv(A);

poly_norm=linear_transform_poly(pdf.poly,Ainv,m(:));
cexp = -log(det(A));
c0 = get_coeff_NDpoly(poly_norm,zeros(1,ndim));
poly_norm = update_or_insert_coeff_NDpoly(poly_norm,zeros(1,ndim),c0+cexp);

normpdf.func=@(x)exp(evaluate_polyND(poly_norm,x));
normpdf.poly=poly_norm;

% now the pdf is normalized
% keyboard
%% dummy MC method
if strcmp(method,'dummyMC')
    Nmc=10000;
    Pint = 0.5*eye(ndim);
    Y = mvnrnd(zeros(1,ndim),Pint,Nmc);
    Pgauss=get_gauss_quad_poly(zeros(1,ndim),Pint);
    
    polly = add_sub_polyND(normpdf.poly,Pgauss,'sub');
    
%     divpdf=mvnpdf(Y,zeros(1,ndim),0.1*eye(ndim));
    s=0;
    for i=1:Nmc
        y= Y(i,:)';
        x = Ainv*y+m;
%         s = s+ F(x)*normpdf.func(y)./divpdf(i);
        d = evaluate_polyND(polly,y);
        d=max(min(d,150),-150);
%         if d==70
%             keyboard
%         end
        
        s = s+ F(x)*exp( d );
    end
    s=s/Nmc;
    
end





f = s;
