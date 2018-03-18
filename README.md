# SYDE 362 - goal of MATLAB component

Our goal is to use the `fmincon` MATLAB function to return an optimal vector x that meets our constraints. We represent the constraints as arguments to the function.

```matlab
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
```

In our case, x is an array of 9 elements. The object looks like:

```matlab
keySet = {'rho', 'Wcb_l', 'Wst_l', 'Wsp_l', 'Wl_l', 'Wcb_k', 'Wst_k', 'Wsp_k', 'Wl_k'};
valueSet = [20, 2, 5, 0.5, 0.25, 1.5, 3, 0.3, 0.15]; % example dummy values
```

## fmincon argument definitions

#### 1. fun is the objective function.

```matlab
objective = @(x) (8*Wcb_l) + (12*Wst_l) + (2*Wsp_l) + (4*Wl_l) + (6*Wcb_k) + (8*Wst_k) + (2*Wsp_k) + (4*Wl_l);
```

where:

```matlab
rho = 0; % material density

Wcb_l = rho*(pi*Rcb_l^2*Lcb_l); % laptop stand crossbar mass
Wst_l = rho*(Lst_l*Tst_l*Hst_l); % laptop stand strut mass
Wsp_l = rho*(Lsp_l*Tsp_l*Hsp_l); % laptop stand support mass
Wl_l = rho*(Tl_l*Ll_l*Hl_l); % laptop stand locking bar mass

Wcb_k = rho*(pi*Rcb_k^2*Lcb_k); % keyboard stand crossbar mass
Wst_k = rho*(Lst_k*Tst_k*Hst_k); % keyboard stand strut mass
Wsp_k = rho*(Lsp_k*Tsp_k*Hsp_k); % keyboard stand support mass
Wl_k = rho*(Tl_k*Ll_k*Hl_k); % keyboard stand locking bar mass
```
#### 2. x0 is our initial guess vector. 

```matlab
x0 = [20, 2, 5, 0.5, 0.25, 1.5, 3, 0.3, 0.15];
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
Aeq = [];
Beq = [];
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







