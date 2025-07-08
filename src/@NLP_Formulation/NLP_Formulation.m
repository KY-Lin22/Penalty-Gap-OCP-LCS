classdef NLP_Formulation < handle
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
    %       c(z, p) >= 0

    properties
        reformulation_strategy char {mustBeMember(reformulation_strategy, {...
            'relaxation',...
            'penalty',...
            'smoothing'
            })} = 'penalty';
        relaxation_problem char {mustBeMember(relaxation_problem, {...
            'Scholtes',...
            'Lin_Fukushima'...
            })} = 'Scholtes' 
        penalty_problem char {mustBeMember(penalty_problem, {...
            'gap_based',...
            'complementarity_based'...
            })} = 'gap_based' 
        smoothing_problem char {mustBeMember(smoothing_problem, {...
            'FB',...
            })} = 'FB'
        D_gap_param_a double {mustBeNonnegative} = 0.9; % D gap function parameters: b > a > 0 (a ref value: a = 0.9, b = 1.1)
        D_gap_param_b double {mustBeNonnegative} = 1.1; % Ref: Theoretical and numerical investigation of the D-gap function   
                                                        % for BVI, 1998, Mathematical Programming, C.Kanzow & M. Fukushima
    end
    properties
        z % symbolic variable, includes all the variable to be optimized,
        p % symbolic variable, including all the problem parameter
        J % symbolic function, cost function 
        h % symbolic function, equality constraint  
        c % symbolic function, inequality constraint      
        Dim % struct, problem dimension record
        
        FuncObj % structure, NLP function object  
    end

    %% Constructor method
    methods
        function self = NLP_Formulation(OCP, Option)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            %% specify properties based on Option
            if isfield(Option, 'reformulation_strategy')
                self.reformulation_strategy = Option.reformulation_strategy;
            end
            if isfield(Option, 'relaxation_problem')
                self.relaxation_problem = Option.relaxation_problem;
            end
            if isfield(Option, 'penalty_problem')
                self.penalty_problem = Option.penalty_problem;
            end
            if isfield(Option, 'smoothing_problem')
                self.smoothing_problem = Option.smoothing_problem;
            end
            if isfield(Option, 'D_gap_param_a')
                self.D_gap_param_a = Option.D_gap_param_a;
            end
            if isfield(Option, 'D_gap_param_b')
                self.D_gap_param_b = Option.D_gap_param_b;
            end
            if self.D_gap_param_b <= self.D_gap_param_a
                error('D gap function parameter should satisfy: b > a > 0')
            end
            
            %% discretize OCP into NLP
            switch self.reformulation_strategy
                case 'relaxation'
                    nlp = self.create_relaxation_NLP(OCP);
                case 'penalty'
                    switch self.penalty_problem
                        case 'gap_based'
                            nlp = self.create_penalty_gap_NLP(OCP);
                        case 'complementarity_based'
                            nlp = self.create_penalty_complementarity_NLP(OCP);
                    end
                case 'smoothing'
                    nlp = self.create_smoothing_NLP(OCP);
            end

            % variable and function
            self.z = nlp.z;   
            self.p = nlp.p;
            self.J = nlp.J;    
            self.h = nlp.h; 
            self.c = nlp.c;       
            % dim
            self.Dim = nlp.Dim;  

            %% create NLP function object
            self.FuncObj = self.create_FuncObj(nlp);
            
            %% display NLP information
            disp('*----------------------------------- NLP Information ------------------------------------*')
            disp('1. complementarity constraint reformulation')
             switch self.reformulation_strategy
                case 'relaxation'
                    disp(['relaxation strategy: .................. ', self.relaxation_problem])
                case 'penalty'
                    disp(['penalty strategy: ..................... ', self.penalty_problem])
                    switch self.penalty_problem
                        case 'gap_based'
                            disp(['D gap parameter (a / b): .............. ',...
                                num2str(self.D_gap_param_a), ' / ', num2str(self.D_gap_param_b)])
                        case 'complementarity_based'
                            
                    end
                case 'smoothing'
                    disp(['smoothing strategy: ................... ', self.smoothing_problem])
            end           
            disp('2. Problem Size')
            disp(['number of decision variable (z): ...... ', num2str(self.Dim.z)])
            disp(['number of equality constraint (h): .... ', num2str(self.Dim.h)])
            disp(['number of inequality constraint (c): .. ', num2str(self.Dim.c)])

        end

    end

    %% Other methods
    methods
        % create relaxation NLP
        nlp = create_relaxation_NLP(self, OCP)

        % create penalty gap NLP
        nlp = create_penalty_gap_NLP(self, OCP)

        D_gap_func = create_D_gap_func(self)

        regular_func = create_regular_func(self)

        [group_func, decouple_func_lambda, decouple_func_eta] ...
            = create_element_wise_concatenation_func(self, OCP)

        % create penalty complementarity NLP
        nlp = create_penalty_complementarity_NLP(self, OCP)

        % create smoothing NLP
        nlp = create_smoothing_NLP(self, OCP)

        % create function object
        FuncObj = create_FuncObj(self, nlp)

    end

end