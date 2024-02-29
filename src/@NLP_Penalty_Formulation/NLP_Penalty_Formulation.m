classdef NLP_Penalty_Formulation < handle
    % formulate a NLP based on the given OCP_LCS and formulation option
    %
    % OCP_LCS has the form:
    %  min  int_0^T 0.5 * (x^T Q_x x + u^T Q_u u + lambda^T Q_lambda lambda) dt,
    %  s.t. Dot{x} = A x + B u + E lambda
    %       eta = C x + D u + F lambda
    %       0 <= lambda \perp eta >= 0 
    %
    % NLP has the form:
    %  min  J(z, p),
    %  s.t. h(z, p) = 0,

    properties
        CHKS_param double {mustBeNonnegative} = 0.001 % used in CHKS smoothing function for max(0, x)
        Huber_param double = 0.1 % used in pseudo Huber loss function, ref: SCP survey 2021, F.Messerer et.al.
        D_gap_param_a double {mustBeNonnegative} = 0.8; % D gap function parameters: b > a > 0 (a ref value: a = 0.9, b = 1.1)
        D_gap_param_b double {mustBeNonnegative} = 1.2; % Ref: Theoretical and numerical investigation of the D-gap function   
                                                        % for BVI, 1998, Mathematical Programming, C.Kanzow & M. Fukushima
    end
    properties
        D_gap_func % function object, D gap function (stage wise)
        Huber_func % function object, pseudo Huber loss function: N -> 1
        Huber_hessian % function object, pseudo Huber loss function hessian
    end
    properties
        z % symbolic variable, includes all the variable to be optimized,
        p % symbolic variable, including all the problem parameter
        J % symbolic function, cost function 
        h % symbolic function, equality constraint    

        J_ocp_hessian % constant matrix, ocp cost hessian
        h_grad % constant matrix, constraint Jacobian

        Dim % struct, problem dimension record
        
        FuncObj % structure, NLP function object  
    end

    %% Constructor method
    methods
        function self = NLP_Penalty_Formulation(OCP, Option)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            %% specify properties based on Option
            if isfield(Option, 'CHKS_param')
                self.CHKS_param = Option.CHKS_param;
            end
            if isfield(Option, 'Huber_param')
                self.Huber_param = Option.Huber_param;
            end
            if isfield(Option, 'D_gap_param_a')
                self.D_gap_param_a = Option.D_gap_param_a;
            end
            if isfield(Option, 'D_gap_param_b')
                self.D_gap_param_b = Option.D_gap_param_b;
            end

            %% specify properties about function object used in NLP reformulation (stage-wise)
            self.D_gap_func = self.create_D_gap_func(OCP);
            [Huber_func, Huber_hessian] = self.create_Huber_func_hessian(OCP);
            self.Huber_func = Huber_func;
            self.Huber_hessian = Huber_hessian;

            %% discretize OCP into NLP
            nlp = self.create_penalty_gap_based_NLP(OCP);
            % variable and function
            self.z = nlp.z;   
            self.p = nlp.p;
            self.J = nlp.J;    
            self.h = nlp.h;   
            % constant jacobian and hessian
            self.J_ocp_hessian = nlp.J_ocp_hessian;
            self.h_grad = nlp.h_grad;
            % dim
            self.Dim = nlp.Dim;  

            %% create NLP function object
            self.FuncObj = self.create_FuncObj(nlp);
            
            %% display NLP information
            disp('*----------------------------------- NLP Information ------------------------------------*')
            disp('1. Problem Parameter')
            disp(['CHKS parameter: ....................... ', num2str(self.CHKS_param, '%10.3e')])
            disp(['Huber parameter: ...................... ', num2str(self.Huber_param, '%10.3e')])
            disp(['D gap parameter (a / b): .............. ', num2str(self.D_gap_param_a), ' / ', num2str(self.D_gap_param_b)])
            disp('2. Problem Size')
            disp(['number of decision variable (z): ...... ', num2str(self.Dim.z)])
            disp(['number of equality constraint (h): .... ', num2str(self.Dim.h)])

        end

    end

    %% Other methods
    methods
        D_gap_func = create_D_gap_func(self, OCP)

        [Huber_func, Huber_hessian] = create_Huber_func_hessian(self, OCP)

        nlp = create_penalty_gap_based_NLP(self, OCP)

        FuncObj = create_FuncObj(self, nlp)

    end

end