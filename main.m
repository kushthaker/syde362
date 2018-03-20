% Objective Function
% FORM: objective = @(x) rho*((8*(pi*Rcb_l^2*Lcb_l)) + (12*(Lst_l*Tst_l*Hst_l)) + (2*(Lsp_l*Tsp_l*Hsp_l)) + (4*Tl_l*Ll_l*Hl_l) + (6*pi*Rcb_k^2*Lcb_k) + (8*(Lst_k*Tst_k*Hst_k)) + (2*(Lsp_k*Tsp_k*Hsp_k)) + (4*Tl_k*Ll_k*Hl_k));

rho = 350;
objective = @(x) x(1)*((8*(pi*x(2)^2*x(3))) + (12*(x(4)*x(5)*x(6))) + (2*(x(7)*x(8)*x(9))) + (4*x(10)*x(11)*x(12)) + (6*pi*x(13)^2*x(14)) + (8*(x(15)*x(16)*x(17))) + (2*(x(18)*x(19)*x(20))) + (4*x(21)*x(22)*x(23)));

% Initial Guess

x0 = [350, 0.0254 / 4, 0.3, 0.404, 0.01, 0.075, 0.404, 0.01, 0.05, 0.3, 0.01, 0.05, 0.0254 / 4, 0.3, 0.404, 0.01, 0.075, 0.404, 0.01, 0.05, 0.3, 0.01, 0.05];

disp(['Initial Guess: ' num2str(objective(x0))]);

% Linear equalities and inequalities

A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper Bounds

lb = zeros([1 23]);
ub = [];

% Non-linear constraints

nonlincon = @nlcon_keyboard; % Both keyboard and laptop constraints need to be in a single function

[x, fval, ef, output]  = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlincon);

% [c, ceq] = nlcon(x);


% keySet = {'rho', 'Rcb_l', 'Lcb_l', 'Lst_l', 'Tst_l', 'Hst_l', 'Lsp_l', 'Tsp_l', 'Hsp_l', 'Tl_l', 'Ll_l', 'Hl_l', 'Rcb_k', 'Lcb_k', 'Lst_k', 'Tst_k', 'Hst_k', 'Lsp_k', 'Tsp_k', 'Hsp_k', 'Tl_k', 'Ll_k', 'Hl_k'};
% valueSet = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]; % dummy values
% mapObj = containers.Map(keySet,valueSet);
% [x, fval, ef, output]  = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlincon);
% disp(x);

disp(['Optimization Result: ' num2str(objective(x))]);