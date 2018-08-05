function margprobs = get_marginalized_exppdf(pdf,mX,PX,margdims,method)
% margdims are the dimensions to be marginalized or removed


margprobs = zeros(size(Xp,1),1);


% keyboard

fixedcols = 1:length(mX);
fixedcols(margdims)=[];

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