function z_Init = polish_initial_guess(self, z_Init)
%UNTITLED27 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

Z_Init = reshape(z_Init, self.NLP.Dim.z_Node(end), self.OCP.nStages);

X_Init      = Z_Init(1 : self.NLP.Dim.z_Node(1), :);
U_Init      = Z_Init(self.NLP.Dim.z_Node(1) + 1 : self.NLP.Dim.z_Node(2), :);
LAMBDA_Init = Z_Init(self.NLP.Dim.z_Node(2) + 1 : self.NLP.Dim.z_Node(3), :);
% ETA_Init    = Z_Init(self.NLP.Dim.z_Node(3) + 1 : self.NLP.Dim.z_Node(4), :);

switch self.Option.polish_initial_guess_method
    case 'eta_g'
        g_FuncObj_map = self.OCP.FuncObj.g.map(self.OCP.nStages);
        g_Init = full(g_FuncObj_map(X_Init, U_Init, LAMBDA_Init));
        ETA_Init = g_Init;
    case 'lambda_posi_eta_zero'
        LAMBDA_Init = abs(LAMBDA_Init);
        ETA_Init = zeros(self.OCP.Dim.lambda, self.OCP.nStages);
end

z_Init = reshape([X_Init; U_Init; LAMBDA_Init; ETA_Init],...
    (self.OCP.Dim.x + self.OCP.Dim.u + self.OCP.Dim.lambda + self.OCP.Dim.lambda) * self.OCP.nStages, 1);
end