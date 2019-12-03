function solveforcoeff_reg(fun,A,b,Aineq,bineq,c0)
% Solve Ac=b with Aineq*c<= bineq
% g is gradient
% J is hessian
% onecol: col-number in A that has all 1s 


if isIllConditioned(A)
    keyboard
    
    Am = zeros(size(A));
    bm = zeros(size(b));
    
    
    for i=1:size(A)
        mx = max(abs(A(i,:)));
        mn = min(abs(A(i,:)));
        mm = mean([mx,mn]);
        
        Am(i,:) = A(i,:)/mm;
        bm(i) = b(i)/mm;
    end
    
    Aineqm = zeros(size(Aineq));
    bineqm = zeros(size(bineq));
    for i=1:size(Aineq)
        mx = max(abs(Aineq(i,:)));
        mn = min(abs(Aineq(i,:)));
        mm = mean([mx,mn]);
        
        Aineqm(i,:) = Aineq(i,:)/mm;
        bineqm(i) = bineq(i)/mm;
    end
    
    
    A = Am;
    b = bm;
    Aineq = Aineqm;
    bineq = bineqm;
end

% options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective');
% options = optimoptions('fminunc','Algorithm','trust-region',...
%     'SpecifyObjectiveGradient',true,'HessianFcn','objective');

% options = optimoptions('fminunc','Algorithm','trust-region','Display','iter');
options = optimoptions('fminunc','MaxIterations',1e3,'Display','iter');

lb= [];
ub= [];
nonlcon = [];
c = fmincon(@(c)norm(c,1),c0,Aineq,bineq,A,b,lb,ub,nonlcon,options);


end
% 
% function [f, g, H] = rosenboth(x)
% % Calculate objective f
% f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
% 
% if nargout > 1 % gradient required
%     g = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
%         200*(x(2)-x(1)^2)];
%     
%     if nargout > 2 % Hessian required
%         H = [1200*x(1)^2-400*x(2)+2, -400*x(1);
%             -400*x(1), 200];  
%     end
% 
% end