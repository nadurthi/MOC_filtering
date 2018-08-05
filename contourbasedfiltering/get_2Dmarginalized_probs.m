function margprobs = get_2Dmarginalized_probs(Xp,n1,n2,X,probs,mX,PX,fullpdf,method)
% Xp is points at which we need probability
% X.probs are the points and weights
% n1,n2 are the dimensions that should remain after marginalization
% fullpdf is corresponding pdf



[xr,xc]=size(Xp);
if xr==1 || xc==1
    Xp=Xp(:)';
end

margprobs = zeros(size(Xp,1),1);


% keyboard

fixedcols = [n1,n2];
margcols = 1:length(mX);
margcols(fixedcols)=[];

if strcmp(method,'dummyMC')
    for i=1:size(Xp,1)
        xfix=Xp(i,:);
        Pnew=get_partial_polyND(fullpdf.poly,xfix,fixedcols);
        pdf.poly = Pnew;
        pdf.func = @(x)evaluate_polyND(Pnew,x);
        
        m=mX(margcols);
        P=PX(margcols,margcols);
        
        
        g = [margcols,fixedcols];
        mmorph = mX(g);
        Pmorph = PX(g,g);
        
        nm = length(margcols);
        m_marg = mmorph(1:nm);
        P_marg = Pmorph(1:nm,1:nm);
        m_fix = mmorph(nm+1:end);
        P_fix = Pmorph(nm+1:end,nm+1:end);
        P_cross = Pmorph(1:nm,nm+1:end);
        
        m= m_marg(:)-P_cross*inv(P_fix)*(xfix(:)-m_fix(:));
        P= P_marg - P_cross*inv(P_fix)*P_cross';
        
%         keyboard
        margprobs(i) = integrate_func_exppdf(@(x)1,pdf,m,P,method);
        
    end
end






end