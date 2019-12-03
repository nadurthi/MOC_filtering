function [X_new, w_new] = pfmeasupdate_custommoc(X_old, w_old, model,Tk1,Tk1step,dttk1, pf, ymeas)
%
% Particle Filter Measurement Update
%
% X_new   - the new set of particles after reweighting and resampling if necessary
%
% Gabriel Terejanu (terejanu@buffalo.edu)

%% ------------------------------------------------------------------------
% weights update
%--------------------------------------------------------------------------
mX_pf = zeros(model.hn, pf.no_particles);
w_new = zeros(size(w_old));
% keyboard
for j = 1 : pf.no_particles
    % get the measurement
%     keyboard
%     mX_pf(:,j) = feval(model.hx, X_old(:,j),model.para_dt);
    mX_pf(:,j) = feval(model.h, X_old(:,j));
    w_new(j) = w_old(j) * getLikelihood(ymeas - mX_pf(:,j), model.R);
end
w_new = w_new./sum(w_new);
X_new = X_old;

[m,P]=MeanCov(X_new',w_new(:));

%% ------------------------------------------------------------------------
% resampling
%--------------------------------------------------------------------------
if (pf.resample) 
	Neff = 1/sum(w_new.^2);
    if (Neff < pf.neff * pf.no_particles)
        I = myresampling(w_new);
        I = round(I);
        for j = 1 : pf.no_particles
            X_new(:,j) = X_old(:,I(j));
        end;    
        % reset weights
        w_new = ones(1,pf.no_particles)/pf.no_particles;
    end
end

%% ------------------------------------------------------------------------
% resampling
%--------------------------------------------------------------------------
if (pf.regularizepf) 
    D=sqrtm(P);
    c=1.1;
    nx=model.fn;
    cnx=pi^(nx/2)/gamma(nx/2+1);
    A=((8/cnx)*(nx+4)*(2*sqrt(pi))^nx)^(1/(nx+4));
    hopt=A*(pf.no_particles)^(1/(nx+4));
    
    
    for j = 1 : pf.no_particles
        while(1)
           g=mvnrnd(zeros(model.fn,1),eye(model.fn)) ;
           u=rand;
           l=EpanechnikovKernel(g,nx,cnx)/(c*mvnpdf(g,zeros(model.fn,1),eye(model.fn)));
           if u<=l
               break
           end
        end
        X_new(:,j) = X_new(:,j) + hopt*D*g;
    end
    
end

end

function K=EpanechnikovKernel(x,nx,cnx)
[r,c]=size(x);
if r==1 || c==1
   x=x(:)'; 
end
[N,nx]=size(x);
% cnx=pi^(nx/2)/gamma(nx/2+1);
K=zeros(N,1);
p=(nx+2)/(2*cnx);
for i=1:N
    if norm(x)<1
        K(i)=p*(1-norm(x(i,:))^2);
    else
        K(i)=0;
    end
end

end
