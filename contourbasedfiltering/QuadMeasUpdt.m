function [xu,Pu]=QuadMeasUpdt(xf,Pf,zm,dt,Tk,model,method)

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
        qd_pts=@(m,P)GH_pts(m,P,4);
    otherwise
        error('smthg is wrong: DONT ask me what')
end

[X,w]=qd_pts(xf,Pf);

Npts = size(X,1);
Z=zeros(Npts,model.hn);
for i=1:Npts
    Z(i,:) = model.h(X(i,:)');
end
[mz,Pz] = MeanCov(Z,w);
Pz = Pz + model.R;

Pcc=CrossCov(X,xf,Z,mz,w);

K=Pcc/Pz;

xu=xf+K*(zm-mz);
Pu=Pf-K*Pz*K';



