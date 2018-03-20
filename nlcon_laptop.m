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

	Wcb_l = rho*(pi*Rcb_l^2*Lcb_l);
    Wst_l = rho*(Lst_l*Tst_l*Hst_l);
    Wsp_l = rho*(Lsp_l*Tsp_l*Hsp_l);
    Wl_l = rho*(Tl_l*Ll_l*Hl_l);


    Ist = (Wst_k^3 * Tst_k) / 12;
	Icb = (pi * Rcb_k^4) / 4;

    dmax = 0;
    dmax_cb = 0;
    Fna = 0;
    theta = 0;
    g = 9.81;

    dmax = (dl_l*(3*(Lcb_k)^3 - 4*(dl_l)^2)*sqrt(8^2 + 8^2)) / 24*Ecb*Icb < 1*e - 3

    Wl_l*0.25 - Wcb_l + 2*(3*Wst_l+Wcb_l)*(1 - 1/(Lst_l*cos(theta)))

    % deflection on top bar - laptop stand

    dmax_cb = g*Wcb_l*dl_l*(3*Lcb_l^2 - 4*dl_l^2) / (4*Ecb*Wl_l*Lcb_l^2) < 1*e - 3

    % beam deflection on strut

    g*Wst_l*(cos(theta) + sin(theta))*Lst_l / (4*Est*Tst_l*Hst_l^3) < 5*e - 3

    % shear force on s pin

    Sz = Wst_l/2 - Wcb_l - Wst_l + dsp_l*Wst_l/(2*Lst_l)
    Pz = Wst_l/2
    Sx = 1*(Sz - Pz)/(tan(theta))

    tau_st_l = sqrt(Sx^2 + Sz^2) / (pi*Rcb_l) < sigma_l / 2

    % Shear force on pin V in lockbar

    Gz = Fna + Wcb_l/2
    Gx = -Sx
    Uz = Gz - Wl_l
    Ux = Gx

    tau_l_l = sqrt(Uz^2 + Ux^2) / pi*r_l^2 < sigma_l / k


















