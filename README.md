# SYDE 362 - Goal of MATLAB component

Our goal is to use the `fmincon` MATLAB function to return an optimal vector x that meets our constraints. We represent the constraints as arguments to the function.

```matlab
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
```

[Read the documentation](https://www.mathworks.com/help/optim/ug/fmincon.html)

In our case, x is an array of 23 elements. The object looks like:

```matlab
% keySet = {'rho', 'Rcb_l', 'Lcb_l', 'Lst_l', 'Tst_l', 'Hst_l', 'Lsp_l', 'Tsp_l', 'Hsp_l', 'Tl_l', 'Ll_l', 'Hl_l', 'Rcb_k', 'Lcb_k', 'Lst_k', 'Tst_k', 'Hst_k', 'Lsp_k', 'Tsp_k', 'Hsp_k', 'Tl_k', 'Ll_k', 'Hl_k'};
% valueSet = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]; % dummy values
```

## naming conventions

The following parameters make up the x vector.

```matlab

 rho = x(1);  % Material Density
 Rcb_l = x(2);  % Radius, Crossbar, Laptop Stand
 Lcb_l = x(3);  % Length, Crossbar, Laptop Stand
 Lst_l = x(4);  % Length, Strut, Laptop Stand
 Tst_l = x(5);  % Thickness, Strut, Laptop Stand
 Hst_l = x(6);  % Height, Strut, Laptop Stand
 Lsp_l = x(7);  % Length, Support, Laptop Stand
 Tsp_l = x(8);  % Thickness, Support, Laptop Stand
 Hsp_l = x(9);  % Height, Support, Laptop Stand
 Tl_l = x(10);  % Thickness, Locking Bar, Laptop Stand
 Ll_l = x(11);  % Length, Locking Bar, Laptop Stand
 Hl_l = x(12);  % Height, Locking Bar, Laptop Stand
 Rcb_k = x(2);  % Radius, Crossbar, Keyboard Stand
 Lcb_k = x(3);  % Length, Crossbar, Keyboard Stand
 Lst_k = x(4);  % Length, Strut, Keyboard Stand
 Tst_k = x(5);  % Thickness, Strut, Keyboard Stand
 Hst_k = x(6);  % Height, Strut, Keyboard Stand
 Lsp_k = x(7);  % Length, Support, Keyboard Stand
 Tsp_k = x(8);  % Thickness, Support, Keyboard Stand
 Hsp_k = x(9);  % Height, Support, Keyboard Stand
 Tl_k = x(10);  % Thickness, Locking Bar, Keyboard Stand
 Ll_k = x(11);  % Length, Locking Bar, Keyboard Stand
 Hl_k = x(12);  % Height, Locking Bar, Keyboard Stand
             
```

The following parameters are used throughout the constraint equations.

```matlab

Wcb_l = rho*(pi*Rcb_l^2*Lcb_l); % laptop stand crossbar mass
Wst_l = rho*(Lst_l*Tst_l*Hst_l); % laptop stand strut mass
Wsp_l = rho*(Lsp_l*Tsp_l*Hsp_l); % laptop stand support mass
Wl_l = rho*(Tl_l*Ll_l*Hl_l); % laptop stand locking bar mass

Wcb_k = rho*(pi*Rcb_k^2*Lcb_k); % keyboard stand crossbar mass
Wst_k = rho*(Lst_k*Tst_k*Hst_k); % keyboard stand strut mass
Wsp_k = rho*(Lsp_k*Tsp_k*Hsp_k); % keyboard stand support mass
Wl_k = rho*(Tl_k*Ll_k*Hl_k); % keyboard stand locking bar mass

```



## fmincon argument definitions 

`x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)`

#### 1. fun is the objective function.

```matlab

% FORM: objective = @(x) rho*((8*(pi*Rcb_l^2*Lcb_l)) + (12*(Lst_l*Tst_l*Hst_l)) + (2*(Lsp_l*Tsp_l*Hsp_l)) + (4*Tl_l*Ll_l*Hl_l) + (6*pi*Rcb_k^2*Lcb_k) + (8*(Lst_k*Tst_k*Hst_k)) + (2*(Lsp_k*Tsp_k*Hsp_k)) + (4*Tl_k*Ll_k*Hl_k));

objective = @(x) x(1)*((8*(pi*x(2)^2*x(3))) + (12*(x(4)*x(5)*x(6))) + (2*(x(7)*x(8)*x(9))) + (4*x(10)*x(11)*x(12)) + (6*pi*x(13)^2*x(14)) + (8*(x(15)*x(16*x(17))) + (2*(x(18)*x(19)*x(20))) + (4*x(21)*x(22)*x(23))));

```

#### 2. x0 is our initial guess vector. 

```matlab
x0 = [20.00, 0.5, 5.00, 3.5, 1.00, 1.02, 1.01, 3.09, 2.55, 0.53, 1.55, 1.23, 2.34, 3.43, 4.23, 1.23, 3.11, 2.23, 2.76, 2.2, 1.17, 1.8, 1.12];
```

#### 3. A,b are matrices representing linear inequality constraints on x

```matlab
A = [];
b = [];
```

#### 4. Aeq,beq are matrices representing linear equality constraints on x

```matlab
Aeq = [];
beq = [];
```

#### 5. lb, ub are lower and upper bounds on x

```matlab
lb = [];
ub = [];
```

#### 6. nonlcon are the nonlinear constraints on x.

We define a single function for our nonlinear constraints to pass to finmincon. The function accepts an array x and returns two arrays, c(x) and ceq(x).

c(x) is the array of nonlinear inequality constraints on x. fmincon attempts to satisfy
* c(x) <= 0 for all entries of c.

ceq(x) is the array of nonlinear equality constraints on x. fmincon attempts to satisfy
* ceq(x) = 0 for all entries of ceq.

For example,

```matlab
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@nlcon)

function [c,ceq] = mycon(x)
c = ...     % Compute nonlinear inequalities at x.
ceq = ...   % Compute nonlinear equalities at x.
```







