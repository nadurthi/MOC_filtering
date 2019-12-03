function S=pmf_sampler(pmfvec,N)
% sample the pmfvec, does not have to be ordered, sample of index is
% returned
pmfvec=pmfvec/sum(pmfvec);
[pmfsort,ind]=sort(pmfvec);
% pmfsort
nids=discretesample(pmfsort, N);

S=ind(nids);








