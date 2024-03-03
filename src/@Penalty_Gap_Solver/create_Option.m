function Option = create_Option(self)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

NLP = self.NLP;
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

%% Option for homotopy
Option.Homotopy.kappa_mu_times = 1.2;
Option.Homotopy.VI_nat_res_tol = 1e-2;

end