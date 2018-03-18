% Objective Function
% FORM: objective = @(x) rho*((8*Wcb_l) + (12*Wst_l) + (2*Wsp_l) + (4*Wl_l) + (6*Wcb_k) + (8*Wst_k) + (2*Wsp_k) + (4*Wl_l));

objective = @(x) x(1)*((8*x(2)) + (12*x(3)) + (2*x(4)) + (4*x(5)) + (6*x(6)) + (8*x(7)) + (2*x(8)) + (4*x(9)));

% Initial Guess

x0 = [20, 2, 5, 0.5, 0.25, 1.5, 3, 0.3, 0.15];

disp(['Initial Guess: ' num2str(objective(x0))]);

% Linear equalities and inequalities

A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper Bounds

lb = [];
ub = [];

% Non-linear constraints

nonlincon = @nlcon; % Both keyboard and laptop constraints need to be in a single function

[x, fval, ef, output]  = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlincon);

[c, ceq] = nlcon(x);


% keySet = {'rho', 'Wcb_l', 'Wst_l', 'Wsp_l', 'Wl_l', 'Wcb_k', 'Wst_k', 'Wsp_k', 'Wl_k'};
% valueSet = [20, 2, 5, 0.5, 0.25, 1.5, 3, 0.3, 0.15]; % dummy values
% mapObj = containers.Map(keySet,valueSet);
% [x, fval, ef, output]  = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlincon);
% disp(x);

disp(['Optimization Result: ' num2str(objective(x))]);