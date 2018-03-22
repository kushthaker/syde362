%% SYDE 362 Group 16 W18 Optimization

% Rcb_l = x(1);
% Tst_l = x(2);
% Hst_l = x(3);
% Tsp_l = x(4);
% Hsp_l = x(5);
% Tl_l = x(6);
% Hl_l = x(7);
% Rcb_k = x(8);
% Tst_k = x(9);
% Hst_k = x(10);
% Tsp_k = x(11);
% Hsp_k = x(12);
% Tl_k = x(13);
% Hl_k = x(14);

% Objective Function
% FORM: objective = @(x) rho*((8*(pi*Rcb_l^2*Lcb_l)) + (12*(Lst_l*Tst_l*Hst_l)) + (2*(Lsp_l*Tsp_l*Hsp_l)) + (4*Tl_l*Ll_l*Hl_l) + (6*pi*Rcb_k^2*Lcb_k) + (8*(Lst_k*Tst_k*Hst_k)) + (2*(Lsp_k*Tsp_k*Hsp_k)) + (4*Tl_k*Ll_k*Hl_k));

% Constant parameters
rho = 350;

Lst_l = 0.384;
Lcb_l = 0.35;
Lsp_l = Lst_l;
Ll_l = Lst_l / 2;

Lst_k = 0.318;
Lcb_k = 0.4;
Lsp_k = Lst_k;
Ll_k = Lst_k / 2;

% Objective function
Mcb_l = @(x) rho * pi * (x(1) ^ 2) * Lcb_l;
Mst_l = @(x) rho * Lst_l * x(2) * x(3);
Msp_l = @(x) rho * Lsp_l * x(4) * x(5);
Ml_l = @(x) rho * Ll_l * x(6) * x(7);

Mcb_k = @(x) rho * pi * (x(8) ^ 2) * Lcb_k;
Mst_k = @(x) rho * Lst_k * x(9) * x(10);
Msp_k = @(x) rho * Lsp_k * x(11) * x(12);
Ml_k = @(x) rho * Ll_k * x(13) * x(14);

M_l = @(x) (8 * Mcb_l(x)) + (12 * Mst_l(x)) + (2 * Msp_l(x)) + (4 * Ml_l(x));
M_k = @(x) (6 * Mcb_k(x)) + (8 * Mst_k(x)) + (2 * Msp_k(x)) + (4 * Ml_k(x));
U = @(x) M_l(x) + M_k(x);
% objective = @(x) rho*((8*(pi*x(2)^2*x(3))) + (12*(x(4)*x(5)*x(6))) + (2*(x(7)*x(8)*x(9))) + (4*x(10)*x(11)*x(12)) + (6*pi*x(13)^2*x(14)) + (8*(x(15)*x(16)*x(17))) + (2*(x(18)*x(19)*x(20))) + (4*x(21)*x(22)*x(23)));

% Initial Guess

x0 = [0.0254 / 4, 0.01, 0.075, 0.01, 0.05, 0.01, 0.05, 0.0254 / 4, 0.01, 0.075, 0.01, 0.05, 0.01, 0.05];
% x0 = [0.0053, 0.0056, 0.0602, 0.0068, 0.048, 0.0068, 0.0480, 0.0062, 0.0063, 0.0636, 0.0066, 0.0493, 0.0085, 0.0453];
% x0 = [0.0047 0.0043 0.0552 0.0052 0.0406 0.0052 0.0406 0.006 0.0048 0.0568 0.005 0.0418 0.0072 0.0378];
% x0 = [0.0047 0.0033 0.0489 0.005 0.0428 0.005 0.0428 0.0058 0.004 0.0537 0.0048 0.0445 0.0072 0.039];

disp(['Initial Guess: ' num2str(U(x0))]);

% Linear equalities and inequalities

A = zeros(6, 14);
b = zeros(6, 1);

A(1, 1) = 4;
A(1, 3) = -1;

A(2, 1) = 4;
A(2, 5) = -1;

A(3, 1) = 4;
A(3, 7) = -1;

A(4, 8) = 4;
A(4, 10) = -1;

A(5, 8) = 4;
A(5, 12) = -1;

A(6, 8) = 4;
A(6, 14) = -1;

Aeq = [];
beq = [];

% Lower and Upper Bounds

lb = 0.001 * ones([1 14]);
ub = [];

% Non-linear constraints

options = optimoptions('fmincon');
options.MaxFunctionEvaluations = 30000;
options.MaxIterations = 10000;
options.Algorithm = 'sqp';

nonlincon = @nlcon; % Both keyboard and laptop constraints need to be in a single function
[x, fval, ef, output]  = fmincon(U, x0, A, b, Aeq, beq, lb, ub, nonlincon, options);

% [c, ceq] = nlcon(x);


% keySet = {'rho', 'Rcb_l', 'Lcb_l', 'Lst_l', 'Tst_l', 'Hst_l', 'Lsp_l', 'Tsp_l', 'Hsp_l', 'Tl_l', 'Ll_l', 'Hl_l', 'Rcb_k', 'Lcb_k', 'Lst_k', 'Tst_k', 'Hst_k', 'Lsp_k', 'Tsp_k', 'Hsp_k', 'Tl_k', 'Ll_k', 'Hl_k'};
% valueSet = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]; % dummy values
% mapObj = containers.Map(keySet,valueSet);
% [x, fval, ef, output]  = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlincon);
% disp(x);

disp(['Optimization Result: ' num2str(U(x))]);

x_best = [0.0041 0.0039 0.0494 0.0055 0.0384 0.0055 0.0385 0.0051 0.0042 0.0485 0.0058 0.0418 0.0081 0.0392];