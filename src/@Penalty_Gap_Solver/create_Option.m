function Option = create_Option(self)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Basic Option
Option.printLevel = 2; % 0: print nothing;  
                       % 1: print results
                       % 2: print results and iteration log (should specified recordLevel as 1)

%% Option for tolerance
Option.maxIterNum = 500;

Option.tol.KKT_error_primal = 1e-6;
Option.tol.KKT_error_dual = 1e-4;
Option.tol.KKT_error_total = 1e-6;
Option.tol.dzNorm = 1e-8;

%% Option for hessian
Option.penalty_hessian_regularization = 1; % 0: without regularization
                                           % 1: with regularization

%% Option for merit line search
Option.LineSearch.betaInit = 1; % initial penalty parameter
Option.LineSearch.rho = 0.1; % desired extend for the negativity of merit function directional derivative
Option.LineSearch.stepSize_Min = 1e-4;
Option.LineSearch.stepSize_DecayRate = 0.5;% choose in (0,1)
Option.LineSearch.nu_D = 1e-4;% desired merit function reduction, default 1e-4

Option.LineSearch.using_Hessian = true;
Option.LineSearch.scaling_constraint_violation = true; %

%% Option for homotopy
Option.Homotopy.kappa_mu_times = 1.2;
Option.Homotopy.VI_nat_res_tol = 1e-2;

end