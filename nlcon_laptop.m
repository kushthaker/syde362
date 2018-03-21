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
    Fkz = 24.45 + Wsp_l;

	Icb = (pi * Rcb_l^4) / 4;
    
    % Will this satisfy the overall sum of z forces?
    Fna = Fkz*0.25 - Wcb_l + 2*(3*Wst_l+Wl_l)*(1 - 1/(Lst_l*cos(theta)));

    % Some concerns here - units don't seem to match up, what are x and z?
    % Also, why aren't there any force components?
    dmax = (dl_l * (3*(Lcb_k)^3 - 4*(dl_l)^2) * sqrt((8 * x^2) + (8 * z^2))) / (24*Ecb*Icb);
    
    c(1) = dmax - 0.001;

    % deflection on top bar - why the g term? may be confusion over whether
    % Wcb is a weight or a mass
    dmax_cb = g*Fkz*dl_l*(3*Lcb_l^2 - 4*dl_l^2) / (4*Ecb*Wl_l*Lcb_l^2);
    c(2) = dmax_cb - 0.001;

    % beam deflection on strut - again, why the g? it doesn't seem to be
    % present in Brian's notes here
    dmax_st = g*Fkz*(cos(theta) + sin(theta)) * (Lst_l^3) / (4*Est*Tst_l*Hst_l^3);
    c(3) = dmax_st - 0.003;

    % shear force on s pin

    Sz = (Fkz/2) - Wcb_l - Wst_l + (dsp_l * Fkz / (2*Lst_l));
    Pz = Fkz/2;
    Sx = (Sz/(tan(theta))) - Pz;

    tau_st_l = sqrt(Sx^2 + Sz^2) / (pi*Rcb_l);
    c(4) = tau_st_l - (sigma_l / 2);

    % Shear force on pin V in lockbar

    Gz = Fna + Wcb_l/2;
    Gx = -Sx;
    Uz = Gz - Wl_l;
    Ux = Gx;

    tau_l_l = sqrt(Uz^2 + Ux^2) / pi*r_l^2;
    c(5) = tau_l_l - (sigma_l / k);
