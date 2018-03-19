objective = @(x) x(1)*x(4)*(x(1)+x(2)+x(3))+x(3)

x0 = [1,5,5,1];

disp(['Initial Objective: ' num2str(objective(x0))])

A = [];
b = [];
Aeq = [];
beq = [];

lb = 1.0* ones(4);
ub = 5.0* ones(4);

nonlincon = @nlcon;

[x, fval, ef, output]  = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlincon);

disp(x);



[c, ceq] = nlcon(x);

disp(['Final Objective: ' num2str(objective(x))])

% keySet = {'rho', 'Rcb_l', 'Lcb_l', 'Lst_l', 'Tst_l', 'Hst_l', 'Lsp_l', 'Tsp_l', 'Hsp_l', 'Tl_l', 'Ll_l', 'Hl_l', 'Rcb_k', 'Lcb_k', 'Lst_k', 'Tst_k', 'Hst_k', 'Lsp_k', 'Tsp_k', 'Hsp_k', 'Tl_k', 'Ll_k', 'Hl_k'};
% valueSet = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]; % dummy values
% mapObj = containers.Map(keySet,valueSet);
% 
% x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23];