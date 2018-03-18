function [c, ceq] = nlcon_laptop_strength(x)


	% keySet = {'rho', 'Rcb_l', 'Lcb_l', 'Lst_l', 'Tst_l', 'Hst_l', 'Lsp_l', 'Tsp_l', 'Hsp_l', 'Tl_l', 'Ll_l', 'Hl_l', 'Rcb_k', 'Lcb_k', 'Lst_k', 'Tst_k', 'Hst_k', 'Lsp_k', 'Tsp_k', 'Hsp_k', 'Tl_k', 'Ll_k', 'Hl_k'};

	% rho = x(1);
	% Rcb_l = x(2);
	% Lcb_l = x(3);
	% Lst_l = x(4);
	% Tst_l = x(5);
	% Hst_l = x(6);
	% Lsp_l = x(7);
	% Tsp_l = x(8);
	% Hsp_l = x(9);
	% Tl_l = x(10);
	% Ll_l = x(11);
	% Hl_l = x(12);
	% Rcb_k = x(13);
	% Lcb_k = x(14);
	% Lst_k = x(15);
	% Tst_k = x(16);
	% Hst_k = x(17);
	% Lsp_k = x(18);
	% Tsp_k = x(19);
	% Hsp_k = x(20);
	% Tl_k = x(21);
	% Ll_k = x(22);
	% Hl_k = x(23);

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











