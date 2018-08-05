function pdf=get_interpolated_pdf(X,probs,mquad,Pquad,Nmoms)


[pdf,~] = get_interp_pdf(X,probs,mquad,Pquad,Nmoms);

% pdf = normalize_exp_pdf(pdf,X,mquad,Pquad,'dummyMC');


