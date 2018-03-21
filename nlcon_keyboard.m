function [c, ceq] = nlcon_keyboard(x)

	% Reference

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
    
    % No nonlinear equalitird
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
    
    % Factor of safety
    K = 2;

	% Crossbar Bar Dimensions - Keyboard Stand

	Rcb_k = x(13);
	Lcb_k = x(14);

	% Strut Dimensions - Keyboard Stand

	Lst_k = x(15);
	Tst_k = x(16);
	Hst_k = x(17);

	% Support Dimensions - Keyboard Stand

	Lsp_k = x(18);
	Tsp_k = x(19);
	Hsp_k = x(20);
	dsp_k = 0;

	% Crossbar Bar Dimensions - Keyboard Stand

	Rcb_k = x(13);
	Lcb_k = x(14);

	% Strut Dimensions - Keyboard Stand

	Lst_k = x(15);
	Tst_k = x(16);
	Hst_k = x(17);

	% Support Dimensions - Keyboard Stand

	Lsp_k = x(18);
	Tsp_k = x(19);
	Hsp_k = x(20);
	dsp_k = 0;

	% Locking Bar Dimensions - Keyboard Stand

	Tl_k = x(21);
	Ll_k = x(22);
	Hl_k = x(23);
	dl_k = (Lcb_k / 2) - Tst_k;

	% Input Params - Keyboard Stand

	Wk = 0.5 * 9.8; % kg * m / s^2
	Fkx = 5.2; % N, 25sin12
	Wsp_k = (Lsp_k * Tsp_k * Hsp_k) * rho;
	Fkz = 24.45 + Wsp_k; % N, 25cos12
	W_kz = Wk + Fkz;
	theta_max = 45;
    theta_min = 5;
    
    theta = theta_max;
    
	k = 2; % factor of safety

	% Objective Function Params - Keyboard Stand

    g = 9.81;
	Wcb_k = (pi * Rcb_k^2 * Lcb_k) * rho * g;
	Wst_k = (Lst_k * Tst_k * Hst_k) * rho * g;
	Wl_k = (Tl_k * Ll_k * Hl_k) * rho * g;
	Wsp_k = (Lsp_k * Tsp_k * Hsp_k) * rho * g;

	% External Reaction Forces - Keyboard Stand

	Ffx = 0.5 * (Fkx);
	Fng = 0.5*(W_kz) + 2*Wst_k + 1.5*Wcb_k + Wl_k + 2*Fkx*tan(theta);
	Fnf = 0.5*(W_kz) + 2*Wst_k + 1.5*Wcb_k + Wl_k - 2*Fkx*tan(theta);

	% Top Crossbar Reaction Forces - Keyboard Stand

	Rcb_az_k = W_kz * (Lsp_k / (2 * Lst_k * cos(theta)));
	Rcb_bz_k = W_kz * (1 - (Lsp_k / (2 * Lst_k * cos(theta))));
	Rcb_bx_k = -Fkx;
    
    Acb_b_k = Tsp_k * (2 * Rcb_k);
    c(7) = (sqrt(Rcb_bx_k^2 + Rcb_bz_k^2) / Acb_b_k) - (sigma_b / K);

	% Top Crossbar Bending Deflection - Keyboard Stand

	Ist = (Wst_k^3 * Tst_k) / 12;
	Icb = (pi * Rcb_k^4) / 4;

	Jcb_az_k = ((Wcb_k * Lcb_k^3) + (2 * Rcb_az_k * ((Lcb_k/2) - dsp_k) * ((3 * Lcb_k^2) - (4 * ((Lcb_k/2) - dsp_k)^2))))/ (48 * Ecb * Icb);
	Jcb_bz_k = ((Wcb_k * Lcb_k^3) + (2 * Rcb_bz_k * ((Lcb_k/2) - dsp_k) * ((3 * Lcb_k^2) - (4 * ((Lcb_k/2) - dsp_k)^2))))/ (48 * Ecb * Icb);
	Jcb_bx_k = (2 * Rcb_bx_k * ((Lcb_k/2) - dsp_k) * ((3 * Lcb_k^2) - (4 * ((Lcb_k/2) - dsp_k)^2))) / (48 * Ecb * Icb);

    % Top crossbar bending deflection
	c(1) = sqrt(Jcb_bz_k^2 + Jcb_bx_k^2) - (0.001 / K);

	% Strut Deflection - Keyboard Stand

	Kx = Rcb_bx_k;
	Kz = Rcb_bz_k;
	Lz = Rcb_az_k;

	Jkm = (((Kx*sin(theta) + Kz*cos(theta)) * (Lst_k/2)^3) / (3*Est*Ist)) + ((Wst_k/2) * cos(theta) * (Lst_k/4)^2 - ((3 * Lst_k/2) - (Lst_k/4))) / (6*Est*Ist);
	Jln = ((Lz*cos(theta) * (Lst_k/2)^3) / (3*Est*Ist)) + (((Wst_k/2) * cos(theta) * (Lst_k/4)^2 * ((3*(Lst_k/2)) - Lst_k/4)) / (6*Est*Ist));
    
    c(2) = Jkm - (0.005 / K);
    c(3) = Jln - (0.005 / K);

	% Reaction forces between Xs - Keyboard Stand

	Dz = (0.5*Wk) + Wst_k + (0.75*Wcb_k) + (Fkx*tan(theta));
	Cz = (0.5*W_kz) + Wst_k + (0.75*Wcb_k) + (Fkx*tan(theta));
	Cx = (Cz + Lz) / tan(theta);
	Dx = -Fkx - Cx;

	% Bottom Strut Reaction Forces - Keyboard Stand

	Az = Wl_k + (0.5*Wcb_k);
	Bz = Az;
	Ez = Bz + Cz - Fng;
	Ex = ((1/tan(theta)) * ((2*Cz) - Ez)) + (2*Cx) - (Hst_k * Ffx / (Lst_k * sin(theta)));

	% Shear Failure - Keyboard Stand

	c(4) = (sqrt(Ex^2 + Ez^2) / pi*r_b^2) - (sigma_b / K);

	% Bottom Crossbar Reaction Forces - Keyboard Stand

	Ax = Ex - Cx - Ffx;
	Bx = Ax;
	Qx = Ax;

	Qz = Az - (0.5 * Wcb_k);
    
    c(8) = (sqrt(Qx^2 + Qz^2) / (2 * Rcb_k * Tl_k)) - (sigma_b / K);

	% Locking Pin Shear Failure - Keyboard Stand

	Ux = Qx;
	Uz = Qz - Wl_k;
	sigma_u = sqrt(Ux^2 + Uz^2) / (pi * r_l^2); 
    c(5) = sigma_u - (sigma_l / k);

    c(9) = (sqrt(Ux^2 + Uz^2) / (2 * Rl_k * Tl_k)) - (sigma_l / K);
    
	% Bottom Crossbar Bending Deflection - Keyboard Stand

	Uaf_z = ((Wcb_k * Lcb_k^3) + (2 * Qz * ((Lcb_k/2) - dl_k) * ((3*Lcb_k^2) - (4*((Lcb_k/2) - dl_k)^2)))) / (48*Ecb*Icb);
	Uaf_x = (2 * Qx * ((Lcb_k/2) - dl_k) * ((3*Lcb_k^2) - (4 * ((Lcb_k/2) - dl_k)^2))) / (48*Ecb*Icb);

	c(6) = sqrt(Uaf_z^2 + Uaf_x^2) - (0.001 / K);















