function normpdf = normalize_exp_pdf(pdf,X,probs,mquad,Pquad,GMM,method)
% X are the points that are tracked

dim=size(X,2);

%estimate normalizing constant
c=1;
%%  method 1 using MC

if strcmp(method,'dummyMC2')
    
    m=mquad;
    Pcov=Pquad;
    
    c=integrate_func_exppdf(@(x)1,pdf,m,Pcov,method);
    disp(['integration constant : ',num2str(c)])
    %     keyboard
end

%% method 2 using mixture of gaussians

if strcmp(method,'GMM_MC')
    if isempty(GMM)
        GMMfitter = GMMFitDataSet(X,probs);
        GMM = GMMfitter.FitGMM_kmeans_optimwt(3);
    end
    C=[];
    prvstd=100;
    for i=1:20
        Nmc = 1000*i;
        Xmc1 = random(MyGmm2MatlabGMM(GMM),Nmc);
        c = mean((pdf.func(Xmc1))./evalGMM(Xmc1,GMM));
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

if strcmp(pdf.pdftype , 'ExpPdf')
    P=pdf.poly;
    cexp = log(1/c);
    c0 = get_coeff_NDpoly(P,zeros(1,dim));
    P = update_or_insert_coeff_NDpoly(P,zeros(1,dim),c0+cexp);

    normpdf = pdf;
    normpdf.poly=P;
    normpdf.func=@(x)exp(evaluate_polyND(P,x));

end

if strcmp(pdf.pdftype , 'HybridPdf')


    normpdf = pdf;
    normpdf.func = @(x)(1/c)*pdf.func(x);

end





end