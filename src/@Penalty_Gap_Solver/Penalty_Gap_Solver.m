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
        KKT_matrix_constant % struct, constant part of the KKT matrix
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
            % properties: KKT_matrix_constant
            J_ocp_hessian = NLP.FuncObj.J_ocp_hessian;
            h_grad = NLP.FuncObj.h_grad;
            
            [i_J_ocp_hessian, j_J_ocp_hessian, s_J_ocp_hessian] = find(J_ocp_hessian);

            [i_h_grad, j_h_grad, s_h_grad] = find(h_grad);
            i_h_grad = i_h_grad + NLP.Dim.z;

            i_h_grad_T = j_h_grad;
            j_h_grad_T = i_h_grad;
            s_h_grad_T = s_h_grad;
            
            i_KKT = [i_J_ocp_hessian; i_h_grad; i_h_grad_T];
            j_KKT = [j_J_ocp_hessian; j_h_grad; j_h_grad_T];
            s_KKT = [s_J_ocp_hessian; s_h_grad; s_h_grad_T];
            KKT_matrix_constant = sparse(i_KKT, j_KKT, s_KKT, NLP.Dim.z + NLP.Dim.h, NLP.Dim.z + NLP.Dim.h, length(s_KKT));

            self.KKT_matrix_constant = KKT_matrix_constant;
            % properties: solver option
            self.Option = self.create_Option();    

            disp('Done!')
        end

    end

    %% Other method
    methods
        % create solver option
        Option = create_Option(self) 

        % solving only single NLP with given p
        [z_Opt, Info] = solve_NLP_single(self, z_Init, p)

        % solve a sequence of NLP in a homotopy manner from p_Init to p_End 
        [z_Opt, Info] = solve_NLP(self, z_Init, p_Init, p_End)

        % merit line search 
        [z_k, Info] = line_search_merit(self, beta, z, dz, p, J, h, J_grad)

        % evaluate natural residual
        natRes = evaluate_natural_residual(self, z_Opt)

        % show result (TO DO)
        show_result(self, Info) 

    end

end