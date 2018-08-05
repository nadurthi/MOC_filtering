function pdf=get_exp_pdf_from_poly(P)
pdf = @(x)exp(evaluate_polyND(P,x));
