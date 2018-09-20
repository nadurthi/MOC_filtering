function [xfp,Pfp]=QuadProp(xf,Pf,dt,Tk,model,method)


switch lower(method)
    case 'ckf'
        qd_pts=@cubature_KF_points;
    case 'ut'
        qd_pts=@(m,P)UT_sigmapoints(m,P,2);
    case 'cut4'
        qd_pts=@conjugate_dir_gausspts;
    case 'cut6'
        qd_pts=@conjugate_dir_gausspts_till_6moment_scheme2;
    case 'cut8'
        qd_pts=@conjugate_dir_gausspts_till_8moment;
    case 'gh'
        qd_pts=@(m,P)GH_pts(m,P,para);
    otherwise
        error('smthg is wrong: DONT ask me what')
end

[X,w]=qd_pts(xf,Pf);

Npts = size(X,1);
Y=zeros(size(X));
for i=1:Npts
    Y(i,:) = model.f(dt,Tk,X(i,:)');
end
[xfp,Pfp] = MeanCov(Y,w);
Pfp = Pfp + model.Q;


