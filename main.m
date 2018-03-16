rho = 0;
Rcb_l = 0;
Lcb_l = 0;
Wst_l = 0;
Hst_l = 0;
Lst_l = 0;
Wsp_l = 0;
Hsp_l = 0;
Lsp_l = 0;
Wl_l = 0;
Hl_l = 0;
Ll_l = 0;
Rcb_p = 0;
Lcb_p = 0;
Wst_p = 0;
Hst_p = 0;
Lst_p = 0;
Wsp_p = 0;
Hsp_p = 0;
Lsp_p = 0;
Tl_p = 0;
Hl_p = 0;
Ll_p = 0;
Wk = 0;
Fkx = 0;

objective = @(x) rho*((8*Wcb_l) + (12*Wst_l) + (2*Wsp_l) + (4*Wl_l) + (6*Wcb_p) + (8*Wst_p) + (2*Wsp_p) + (4*Wl_l)); 

keySet = {'rho', 'Wcb_l', 'Wst_l', 'Wsp_l', 'Wl_l', 'Wcb_p', 'Wst_p', 'Wsp_p', 'Wl_p'};
valueSet = [20, 2, 5, 0.5, 0.25, 1.5, 3, 0.3, 0.15];
mapObj = containers.Map(keySet,valueSet);

W_p = Wk + Fkx;

Wcb_l = pi*Rcb_l^2*Lcb_l;
Wst_l = Lst_l*Tst_l*Hst_l;
Wsp_l = Lsp_l*Tsp_l*Hsp_l;
Wl_l = Tl_l*Ll_l*Hl_l;

Wcb_p = pi*Rcb_p^2*Lcb_p;
Wst_p = Lst_p*Tst_p*Hst_p;
Wsp_p = Lsp_p*Tsp_p*Hsp_p;
Wl_p = Tl_p*Ll_p*Hl_p;

Ffx = 0.5*(Fkx);
Fng = 0.5*(W_p) + 2*Wst_p + 1.5*Wcb_p + Wl_p + 2*Fkx*tan(theta);
Fnf = 0.5*(W_p) + 2*Wst_p + 1.5*Wcb_p + Wl_p - 2*Fkx*tan(theta);