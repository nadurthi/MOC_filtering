function dr=getALLcorners(mnb,mxb)
mnb=mnb(:)';
mxb=mxb(:)';
dim = length(mnb);
dr=prod_conjugate_dir(dim);
dr=(dr+1)/2;

N=size(dr,1);
% get all the 4 corners of a box given the principle diagonal corners
dr=dr.*repmat(mxb-mnb,N,1);
dr=dr+repmat(mnb,N,1);

