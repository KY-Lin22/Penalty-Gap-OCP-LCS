function z_Init = polish_initial_guess(self, z_Init)
%UNTITLED27 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

Z_Init = reshape(z_Init, self.NLP.Dim.z_Node(end), self.OCP.nStages);

X_Init      = Z_Init(1 : self.NLP.Dim.z_Node(1), :);
U_Init      = Z_Init(self.NLP.Dim.z_Node(1) + 1 : self.NLP.Dim.z_Node(2), :);
XI_Init     = Z_Init(self.NLP.Dim.z_Node(2) + 1 : self.NLP.Dim.z_Node(4), :);

% decouple Xi
[group_func, decouple_func_lambda, ~] = self.NLP.create_element_wise_concatenation_func(self.OCP);

decouple_func_lambda_map = decouple_func_lambda.map(self.OCP.nStages);
LAMBDA_Init = full(decouple_func_lambda_map(XI_Init));

switch self.Option.polish_initial_guess_method
    case 'eta_g'
        g_FuncObj_map = self.OCP.FuncObj.g.map(self.OCP.nStages);
        g_Init = full(g_FuncObj_map(X_Init, U_Init, LAMBDA_Init));
        ETA_Init = g_Init;
    case 'lambda_posi_eta_zero'
        LAMBDA_Init = abs(LAMBDA_Init);
        ETA_Init = zeros(self.OCP.Dim.lambda, self.OCP.nStages);
end
group_func_map = group_func.map(self.OCP.nStages);
XI_Init = full(group_func_map(LAMBDA_Init, ETA_Init));
z_Init = reshape([X_Init; U_Init; XI_Init],...
    (self.OCP.Dim.x + self.OCP.Dim.u + self.OCP.Dim.lambda + self.OCP.Dim.lambda) * self.OCP.nStages, 1);
end