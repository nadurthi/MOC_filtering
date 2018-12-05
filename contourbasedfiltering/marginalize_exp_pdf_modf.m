function margval = marginalize_exp_pdf_modf(xeval,states2remove,pdfnorm,X,probs,mquad,Pquad,method)
% xeval evaluate at this point
% first remove the states2remove

dim=size(X,2);
states = 1:dim;
states2keep = states;
states2keep(states2remove)=[];
%estimate normalizing constant
c=1;
%%  method 1 using MC

if strcmp(method,'dummyMC2')
    
    m=mquad;
    Pcov=Pquad;
    
    c=integrate_func_exppdf(@(x)1,pdfnorm,m,Pcov,method);
    disp(['integration constant : ',num2str(c)])
    %     keyboard
end

%% method 2 using mixture of gaussians

if strcmp(method,'GMM_MC')
   
    GMMmarg=pdfnorm.GMMHull.GMMhull;
    for i=1:GMMmarg.Ngcomp
        GMMmarg.mx{i} = GMMmarg.mx{i}(states2remove);
        GMMmarg.Px{i} = GMMmarg.Px{i}(states2remove,states2remove);
    end
    
    C=[];
    prvstd=100;
    for i=[5]
        Nmc = 1000*i;
        Xmc1 = random(MyGmm2MatlabGMM(GMMmarg),Nmc);
        XXX = zeros(Nmc,dim);
        XXX(:,states2remove) = Xmc1;
        XXX(:,states2keep) = repmat(xeval(:)',Nmc,1);
        
        c = mean((pdfnorm.func(XXX))./GaussSumMix(Xmc1,GMMmarg));
        C = [C,c];
        sdC = std(C);
        if abs(sdC - prvstd)/prvstd < 0.2
            disp(['normalization integral converged with #samples = ',num2str(Nmc)])
            break
        end
        prvstd = sdC;
    end
    disp(['integration constant : ',num2str(c)])
    %     keyboard
end




%% method 3 do not use X, but do quadratic approximation to compute mean and covariance, and then do sampling integration






%% now see to append the normalizing constant into the exp pdf

margval = c;


end