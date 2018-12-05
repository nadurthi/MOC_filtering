function X = resampling_GMMfut_GH_2D(X,p)
   
    
    X=[];
    for i=1:Ngcomp
       if obj.GMM.w(i)>0.7
           Npara=7;
           c=0.2;
       elseif obj.GMM.w(i)>0.5 && obj.GMM.w(i)<=0.7  
           Npara=7;
           c=0.2;
       elseif obj.GMM.w(i)>0.3 && obj.GMM.w(i)<=0.5
           Npara=7;
           c=0.2;
       elseif obj.GMM.w(i)>0.1 && obj.GMM.w(i)<=0.3
           Npara=5;
           c=0.2;
       elseif obj.GMM.w(i)>1e-2 && obj.GMM.w(i)<=0.1
           Npara=4;
           c=0.2;
       elseif obj.GMM.w(i)>1e-5 && obj.GMM.w(i)<=1e-2
           Npara=3;
           c=0.2;
       elseif obj.GMM.w(i)<1e-5
           Npara=0;
           c=0.1;
       end
       if Npara==0
           x=[];
       else
            x = GH_pts(obj.GMM.mx{i},c^2*obj.GMM.Px{i},Npara);
       end
       X=vertcat(X,x);
    end

end