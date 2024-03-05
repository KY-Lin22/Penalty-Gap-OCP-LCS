function Option = create_Option(self)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Basic Option
Option.printLevel = 2; % 0: print nothing;  
                       % 1: print results
                       % 2: print results and iteration log (should specified recordLevel as 1)

%% Option for tolerance
Option.maxIterNum = 500;

Option.tol.KKT_error_primal = 1e-8;
Option.tol.KKT_error_dual = 1e-6;
Option.tol.KKT_error_total = 1e-8;
Option.tol.dzNorm = 1e-10;

%% Option for initial guess
Option.polish_initial_guess_method = 'lambda_posi_eta_zero'; % 'eta_g', 'lambda_posi_eta_zero'

%% Option for hessian
Option.penalty_hessian_regularization = 1; % 0: without regularization
                                           % 1: with regularization
Option.penalty_hessian_additional_regular_param = 1e-8;

%% Option for QP solver to evaluate search direction
Option.qpSolver = 'direct_sparse'; % 'fbstab_sparse', 'osqp', 'direct_sparse'
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

% initial guess for fbstab_sparse 
Option.fbstab_x0.z = zeros(self.NLP.Dim.z, 1); % primal variable
Option.fbstab_x0.l = zeros(self.NLP.Dim.h, 1); % dual variable for equality constraint h
Option.fbstab_x0.v = zeros(1, 1); % dual variable for inequality constraint c

% option for osqp
osqp_solver = osqp;
Option.osqp_options = osqp_solver.default_settings(); 
Option.osqp_options.verbose = 0;
Option.osqp_options.scaling = 0;
Option.osqp_options.max_iter = 5000;

% initial guess for osqp
Option.osqp_x0 = zeros(self.NLP.Dim.z, 1); % primal variable
Option.osqp_y0 = zeros(self.NLP.Dim.h, 1); % dual variable for equality constraint h and inequality constraint c 

%% Option for merit line search
Option.LineSearch.betaInit = 1; % initial penalty parameter
Option.LineSearch.rho = 0.1; % desired extend for the negativity of merit function directional derivative
Option.LineSearch.stepSize_Min = 1e-4;
Option.LineSearch.stepSize_DecayRate = 0.5;% choose in (0,1)
Option.LineSearch.nu_D = 1e-4;% desired merit function reduction, default 1e-4

Option.LineSearch.scaling_constraint_violation = true; %

%% Option for homotopy
Option.Homotopy.kappa_mu_times = 1.2;
Option.Homotopy.VI_nat_res_tol = 1e-2;

end