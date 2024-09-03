classdef OCP_Formulation < handle
    %create an OCP_LCS with the form of:
    %  min  int_0^T 0.5 * (x^T Q_x x + u^T Q_u u + lambda^T Q_lambda lambda) dt,
    %  s.t. Dot{x} = A x + B u + E lambda
    %       eta = C x + D u + F lambda
    %       0 <= lambda \perp eta >= 0

    properties
        TimeHorizon % time horizon
        nStages % number of discretized stage
        timeStep % discretization time step

        x0 % initial state

        x % differentiable state
        u % control input
        lambda % algebraic variable  

        Q_x % stage cost weighting matrix Q_x
        Q_u % stage cost weighting matrix Q_u
        Q_lambda % stage cost weighting matrix Q_lambda
        
        L_S % stage cost  

        A % ODE r.h.s function matrix A
        B % ODE r.h.s function matrix B
        E % ODE r.h.s function matrix E
        C % complementarity function matrix C
        D % complementarity function matrix D
        F % complementarity function matrix F        
        
        f % ODE r.h.s function
        g % complementarity function
        
        Dim % variable dimemsion record   
        FuncObj % CasADi function object
    end

    %% Constructor method  
    methods
        function self = OCP_Formulation(...
                TimeHorizon, nStages, timeStep,...
                x0, ...
                x, u, lambda,...
                Q_x, Q_u, Q_lambda,...
                A, B, E,...
                C, D, F)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            disp('creating OCP of LCS...')
            %% formulate OCP_LCS
            % time parameter
            self.TimeHorizon = TimeHorizon;
            self.nStages = nStages;
            self.timeStep = timeStep;
            % initial state
            self.x0 = x0;   
            % variable 
            self.x = x;
            self.u = u;
            self.lambda = lambda; 
            % cost function
            self.Q_x = Q_x;
            self.Q_u = Q_u;
            self.Q_lambda = Q_lambda;
            self.L_S = 0.5 * (x' * Q_x * x + u' * Q_u * u + lambda' * Q_lambda * lambda);
            % linear complementarity system
            self.A = A;
            self.B = B;
            self.E = E;
            self.C = C;
            self.D = D;
            self.F = F;
            self.f = A * x + B * u + E * lambda;
            self.g = C * x + D * u + F * lambda;
            % dim record
            self.Dim = struct('x', size(x, 1), 'u', size(u, 1), 'lambda', size(lambda, 1)); 

            % function object
            self.FuncObj = self.create_FuncObj(); 

            %% display OCP_LCS informulation
            disp('*---------------------------------- OCP Information -----------------------------------*')
            disp('1. time parameter')
            disp(['time horizon: ............................... ', num2str(self.TimeHorizon)])
            disp(['discretization stage: ....................... ', num2str(self.nStages)])
            disp(['time step: .................................. ', num2str(self.timeStep)])            
            disp('2. problem size')
            disp(['number of state variable (x): ............... ', num2str(self.Dim.x)])
            disp(['number of control variable (u): ............. ', num2str(self.Dim.u)])
            disp(['number of algebraic variable (lambda): ...... ', num2str(self.Dim.lambda)])         

            disp('Done!')
        end

    end

    %% other method
    methods
        FuncObj = create_FuncObj(self)
    end
end