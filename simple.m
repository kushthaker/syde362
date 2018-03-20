% FORM: objective = @(x) rho*((8*(pi*Rcb_l^2*Lcb_l)) + (12*(Lst_l*Tst_l*Hst_l)) + (2*(Lsp_l*Tsp_l*Hsp_l)) + (4*Tl_l*Ll_l*Hl_l) + (6*pi*Rcb_k^2*Lcb_k) + (8*(Lst_k*Tst_k*Hst_k)) + (2*(Lsp_k*Tsp_k*Hsp_k)) + (4*Tl_k*Ll_k*Hl_k));
objective = @(x) x(1)*((6*pi*x(13)^2*x(14)) + (8*(x(15)*x(16*x(17))) + (2*(x(18)*x(19)*x(20))) + (4*x(21)*x(22)*x(23))));


x0 = [20.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 4.23, 1.23, 3.11, 2.23, 2.76, 2.22, 1.17, 1.82, 1.12];

disp(['Initial Guess: ' num2str(objective(x0))]);

% Linear equalities and inequalities

A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper Bounds

lb = [];
ub = [];