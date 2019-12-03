function [X_new] = pftimeupdate_custommoc(X_old, model,Tk1,Tk1step,dtk1)
%
% Particle Filter Time Update - DISCRETE DYNAMICAL SYSTEM
%
% X_new   - the new set of particles after propagation and noise perturb
%
% Gabriel Terejanu (terejanu@buffalo.edu)

no_particles = size(X_old, 2);
X_new = zeros(size(X_old));
% keyboard
for i = 1 : no_particles
    %     [t yend y] = feval(model.fx, X_old(:,i),time.tspan(k-measTs), time.tspan(k));
    %  yend = feval(model.fx, X_old(:,i),model.para_dt);
    
        yend=model.f(dtk1,Tk1,X_old(:,i));

    [r,c]=size(yend);
    if r==model.fn
        yend=yend';
    end
    X_new(:,i) = yend' + model.sQ * randn(model.fn,1);
    %     Y(i,:,:)=y;
    %     i
end

% keyboard
% for j=1
%     for indt=1:length(t)
%         Momx1(indt,j)=sum(Y(:,indt,1))/no_particles;
%         Momx2(indt,j)=sum(Y(:,indt,2))/no_particles;
%         Momp1(indt,j)=sum(Y(:,indt,3))/no_particles;
%         Momp2(indt,j)=sum(Y(:,indt,4))/no_particles;
%     end
% end
%
% for j=2:Nm
%     for indt=1:length(t)
%             Momx1(indt,j)=sum((Y(:,indt,1)-Momx1(indt,1)).^j)/no_particles;
%             Momx2(indt,j)=sum((Y(:,indt,2)-Momx2(indt,1)).^j)/no_particles;
%             Momp1(indt,j)=sum((Y(:,indt,3)-Momp1(indt,1)).^j)/no_particles;
%             Momp2(indt,j)=sum((Y(:,indt,4)-Momp2(indt,1)).^j)/no_particles;
%     end
%     j
% end
