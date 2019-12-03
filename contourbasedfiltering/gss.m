function [lb,ub]=gss(f, a, b, tol)
invphi = (sqrt(5) - 1) / 2;  %# 1 / phi
invphi2 = (3 - sqrt(5)) / 2;  %# 1 / phi^2

%    Given a function f with a single local minimum in
%    the interval [a,b], gss returns a subset interval
%    [c,d] that contains the minimum with d-c <= tol.


    a = min(a, b);
    b = max(a, b);
    h = b - a;
    if h <= tol
        lb=a;
        ub=b;
        return 
    end
%     # Required steps to achieve tolerance
    n = round(ceil(log(tol / h) / log(invphi)));

    c = a + invphi2 * h;
    d = a + invphi * h;
    yc = f(c);
    yd = f(d);

    for k =1:n-1
        if yc < yd
            b = d;
            d = c;
            yd = yc;
            h = invphi * h;
            c = a + invphi2 * h;
            yc = f(c);
        else
            a = c;
            c = d;
            yc = yd;
            h = invphi * h;
            d = a + invphi * h;
            yd = f(d);
        end
    end
    if yc < yd
        lb=a;
        ub=d;
        return ;
    else
        lb=c;
        ub=b;
        return;
    end