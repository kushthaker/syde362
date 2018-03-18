function [c, ceq] = nlcon_laptop_strength(x)

	% Material Properties

	rho = 0;
	Est = 0;
	Ecb = 0;

	% Bearing and Pin Dimensions

	r_b = 0;
	sigma_b = 0;
	r_l = 0;
	sigma_l = 0;

	% Crossbar Dimensions - Laptop Stand

	Rcb_l = 0;
	Lcb_l = 0;

	% Strut Dimensions - Laptop Stand

	Tst_l = 0;
	Hst_l = 0;
	Lst_l = 0;

	% Support Dimensions - Laptop Stand

	Tsp_l = 0;
	Hsp_l = 0;
	Lsp_l = 0;
	dsp_l = 0;

	% Locking Bar Dimensions - Laptop Stand

	Tl_l = 0;
	Hl_l = 0;
	Ll_l = 0;
	dl_l = 0;

	% Objective Function Params - Laptop Stand

	Wcb_l = (pi*Rcb_l^2*Lcb_l);
    Wst_l = (Lst_l*Tst_l*Hst_l);
    Wsp_l = (Lsp_l*Tsp_l*Hsp_l);
    Wl_l = (Tl_l*Ll_l*Hl_l);











