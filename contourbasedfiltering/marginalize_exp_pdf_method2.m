function margprob = marginalize_exp_pdf_method2(fixedxs,fixedstates,normpdffunc,dim)
% xeval evaluate at this point
% first remove the states2remove
% 
[r,c]=size(fixedxs);
if r==1||c==1
    fixedxs=fixedxs(:)';
end

[N,~]=size(fixedxs);

margstates = setdiff(1:dim,fixedstates);

%%  method 1 using MC
xl = -1*ones(1,length(margstates));
xu =  1*ones(1,length(margstates));
vv = prod(xu-xl);
margprob = zeros(N,1);

for i=1:N
    fx = fixedxs(i,:);
    
    prev=0.0001;
    for Ninteg=9
        [xint,wint] = GLeg_pts(Ninteg*ones(1,length(margstates)), xl, xu);

        y=fixedcoordpdfexp(normpdffunc,fixedstates,xint,fx,dim);
        m=y(:)'*wint(:);
%         [m,prev,abs(m-prev)/prev]
        if abs(m-prev)/prev<1e-1
            break
        end
        prev=m;
    end
    
    margprob(i)=vv*m;



end