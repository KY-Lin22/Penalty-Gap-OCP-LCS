classdef Penalty_Gap_Solver < handle
    %Implementation of penalty gap method which solves the following NLP problem:
    %  
    %  min  J(z, p),
    %  s.t. h(z, p) = 0,
    % where: z is the variable, 
    %        p is the parameter, 
    %        J is the cost, and h is the constraint arranged in a stagewise manner.
    properties
        OCP % struct, optimal control problem
        NLP % struct, nonlinear programming problem (discretized OCP)
        Option % struct, IPOPT solver option
    end
    %% Constructor method  
    methods
        function self = Penalty_Gap_Solver(OCP, NLP)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            disp('creating solver...')
             % properties: OCP and NLP
            self.OCP = OCP;
            self.NLP = NLP;            
            
            % properties: solver option
            self.Option = self.create_Option();    

            disp('Done!')
        end

    end

    %% Other method
    methods
        % create solver option
        Option = create_Option(self) 

        % create initial guess (TODO)
        z_Init = create_initial_guess(self)

        % solving only single NLP with given p
        [z_Opt, Info] = solve_NLP_single(self, z_Init, p)

        % solve a sequence of NLP in a homotopy manner from p_Init to p_End 
        [z_Opt, Info] = solve_NLP(self, z_Init, p_Init, p_End)

        % evaluate search direction
        [dz, gamma_h_k, Info] = evaluate_search_direction(self, h, J_grad, h_grad, LAG_hessian)

        % merit line search 
        [z_k, Info] = line_search_merit(self, beta, z, dz, p, J, h, J_grad, LAG_hessian)

        % evaluate natural residual
        natRes = evaluate_natural_residual(self, z_Opt)

        % show result (TO DO)
        show_result(self, Info) 

    end

end