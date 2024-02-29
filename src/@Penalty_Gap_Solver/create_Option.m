function Option = create_Option(self)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% Basic Option
Option.printLevel = 2; % 0: print nothing;  
                       % 1: print results
                       % 2: print results and iteration log (should specified recordLevel as 1)
Option.recordLevel = 1; % 0: record time 
                        % 1: record time and log

%% Option for tolerance
Option.maxIterNum = 100;

Option.KKT_scaling_max = 100;
Option.tol.KKT_error_primal = 1e-10;
Option.tol.KKT_error_dual = 1e-10;
Option.tol.KKT_error_total = 1e-8;
Option.tol.dzNorm = 1e-8;

%% Option for Hessian regularization



%% Option for QP solver to evaluate search direction
Option.qpSolver = 'fbstab_sparse'; % 'fbstab_sparse', 'osqp'
% option for fbstab_sparse
Option.fbstab_options = fbstab_options();
Option.fbstab_options.display_level = 0;
Option.fbstab_options.max_newton_iters = 1000;
Option.fbstab_options.solver_mode = 3; % fbstab_sparse only:
                                       %   solver_mode: (default: 3, seems to be the most efficient choice)
                                       %     1: solve the asymmetric newton step system with sparse LU
                                       %     2: solve the quasidefinite reduced newton system with LDL'
                                       %     3: solve the normal equations with LL'
Option.fbstab_options.use_ordering = true; % default: true

% option for osqp
osqp_solver = osqp;
Option.osqp_options = osqp_solver.default_settings(); 
Option.osqp_options.verbose = 0;
Option.osqp_options.scaling = 0;
Option.osqp_options.max_iter = 5000;

%% Option for merit line search
Option.LineSearch.betaInit = 1; % initial penalty parameter
Option.LineSearch.rho = 0.1; % desired extend for the negativity of merit function directional derivative
Option.LineSearch.stepSize_Min = 0.001;
Option.LineSearch.stepSize_DecayRate = 0.5;% choose in (0,1)
Option.LineSearch.nu_D = 1e-4;% desired merit function reduction, default 1e-4

Option.LineSearch.scaling_constraint_violation = false; %

%% Option for homotopy
Option.Homotopy.kappa_mu_times = 1.2;

Option.Homotopy.VI_nat_res_tol = 1e-2;

end