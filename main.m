% Material Properties

rho = 0;
Est = 0;
Ecb = 0;

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

% Crossbar Bar Dimensions - Keyboard Stand

Rcb_k = 0;
Lcb_k = 0;

% Strut Dimensions - Keyboard Stand

Tst_k = 0;
Hst_k = 0;
Lst_k = 0;

% Support Dimensions - Keyboard Stand

Tsp_k = 0;
Hsp_k = 0;
Lsp_k = 0;
dsp_k = 0;

% Locking Bar Dimensions - Keyboard Stand

Tl_k = 0;
Hl_k = 0;
Ll_k = 0;

% Input Params

Wk = 0;
Fkx = 0;
Fkz = 0;
W_kz = Wk + Fkz;
theta = 0;

% Objective Function Params

Wcb_l = rho*(pi*Rcb_l^2*Lcb_l);
Wst_l = rho*(Lst_l*Tst_l*Hst_l);
Wsp_l = rho*(Lsp_l*Tsp_l*Hsp_l);
Wl_l = rho*(Tl_l*Ll_l*Hl_l);

Wcb_k = rho*(pi*Rcb_k^2*Lcb_k);
Wst_k = rho*(Lst_k*Tst_k*Hst_k);
Wsp_k = rho*(Lsp_k*Tsp_k*Hsp_k);
Wl_k = rho*(Tl_k*Ll_k*Hl_k);

% Objective Function

objective = @(x) (8*Wcb_l) + (12*Wst_l) + (2*Wsp_l) + (4*Wl_l) + (6*Wcb_k) + (8*Wst_k) + (2*Wsp_k) + (4*Wl_l);

% Initial Guess Map Object

keySet = {'rho', 'Wcb_l', 'Wst_l', 'Wsp_l', 'Wl_l', 'Wcb_k', 'Wst_k', 'Wsp_k', 'Wl_k'};
valueSet = [20, 2, 5, 0.5, 0.25, 1.5, 3, 0.3, 0.15];
mapObj = containers.Map(keySet,valueSet);

% External Reaction Forces

Ffx = 0.5*(Fkx);
Fng = 0.5*(W_kz) + 2*Wst_k + 1.5*Wcb_k + Wl_k + 2*Fkx*tan(theta);
Fnf = 0.5*(W_kz) + 2*Wst_k + 1.5*Wcb_k + Wl_k - 2*Fkx*tan(theta);

% Top Crossbar Reaction Forces

Rcb_az_k = W_kz*(Lsp_k/(Lst_k*cos(theta)));
Rcb_bz_k = W_kz*(1 - (Lsp_k/(Lst_k*cos(theta))));
Rcb_bx_k = -Fkx;

% Top Crossbar Reaction Forces

Ist = (Wst_k^3*Tst_k)/12;
Icb = (pi*Rcb_k)/4;

Jcb_az_k = ((Wcb_k*Lcb_k)^3 + 2*Rcb_az_k*((Lcb_k/2) - dsp_k)*(3*Lcb_k^2 - 4*((Lcb_k/2) - dsp_k)^2))/48*Ecb*Icb;
Jcb_bz_k = ((Wcb_k*Lcb_k)^3 + 2*Rcb_bz_k*((Lcb_k/2) - dsp_k)*(3*Lcb_k^2 - 4*((Lcb_k/2) - dsp_k)^2))/48*Ecb*Icb;
Jcb_bx_k = (2*Rcb_bx_k*((Lcb_k/2) - dsp_k)*(3*Lcb_k^2 - 4*((Lcb_k/2) - dsp_k)^2))/48*Ecb*Icb;
