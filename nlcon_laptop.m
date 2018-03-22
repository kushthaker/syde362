function [c, ceq] = nlcon_laptop(x)


	% keySet = {'rho', 'Rcb_l', 'Lcb_l', 'Lst_l', 'Tst_l', 'Hst_l', 'Lsp_l', 'Tsp_l', 'Hsp_l', 'Tl_l', 'Ll_l', 'Hl_l', 'Rcb_k', 'Lcb_k', 'Lst_k', 'Tst_k', 'Hst_k', 'Lsp_k', 'Tsp_k', 'Hsp_k', 'Tl_k', 'Ll_k', 'Hl_k'};

    % no nonlinear equalities
    ceq = [];
    
	% Material Properties
	rho = 350; % (x(1)) 350 kg / m^3
	Est = 8500 * 1000; % kPa
	Ecb = 8500 * 1000; % kPa

	% Bearing and Pin Dimensions
    r_b = .0254 / 2 / 2; % half inch diameter
	sigma_b = 6200; % kPa
	r_l = .0254 / 2 / 2; % half inch diameter
	sigma_l = 6200; % kPa

	% Crossbar Dimensions - Laptop Stand
    Rcb_l = x(2);
	Lcb_l = x(3);
    
	% Strut Dimensions - Laptop Stand
	Lst_l = x(4);
	Tst_l = x(5);
	Hst_l = x(6);

	% Support Dimensions - Laptop Stand
	Lsp_l = x(7);
	Tsp_l = x(8);
	Hsp_l = x(9);
	dsp_l = 0;

	% Locking Bar Dimensions - Laptop Stand
	Tl_l = x(10);
	Ll_l = x(11);
	Hl_l = x(12);
	dl_l = (Lcb_l / 2) - Tst_l;
    
	% Objective Function Params - Laptop Stand
    g = 9.81;
    
	theta_max = 45;
    theta_min = 5;
    
    theta = theta_max;

	Wcb_l = rho * (pi * Rcb_l^2 * Lcb_l) * g;
    Wst_l = rho * (Lst_l * Tst_l * Hst_l) * g;
    Wsp_l = rho * (Lsp_l * Tsp_l * Hsp_l) * g;
    Wl_l = rho * (Tl_l * Ll_l * Hl_l) * g;
    Flz = 24.45 + Wsp_l;

	Icb = (pi * Rcb_l^4) / 4;
    
    Fna = (Flz / 4) + (2 * Wcb_l) + (3*Wst_l) + Wl_l;
    
    dmax_cb = Flz*dl_l*(3*Lcb_l^2 - 4*dl_l^2) / (48 * Ecb * (pi * Rcb_l^4 / 4));
    c(2) = dmax_cb - (0.001 / K);

    % beam deflection on strut - again, why the g? it doesn't seem to be
    % present in Brian's notes here
    dmax_st = Flz*(cos(theta) + sin(theta)) * (Lst_l^3) / (48 * Est * (Tst_l * Hst_l^3 / 12));
    c(3) = dmax_st - (0.003 / K);

    % shear force on s pin

    Sz = (Flz/2) - Wcb_l - Wst_l + (dsp_l * Flz / (2*Lst_l));
    Pz = Flz/2;
    Sx = (Sz/(tan(theta))) - Pz;

    tau_st_l = sqrt(Sx^2 + Sz^2) / (pi*Rcb_l^2);
    c(4) = tau_st_l - (sigma_l / K);

    % Shear force on pin V in lockbar

    Gz = Fna + Wcb_l/2;
    Gx = -Sx;
    Uz = Gz - Wl_l;
    Ux = Gx;

    tau_l_l = sqrt(Uz^2 + Ux^2) / (pi*r_l^2);
    c(5) = tau_l_l - (sigma_l / K);
    
    dmax = (dl_l * (3*(Lcb_l)^2 - 4*(dl_l)^2) * sqrt((Gx^2) + (Gz^2))) / (24*Ecb*Icb);
    c(1) = dmax - (0.001 / K);
