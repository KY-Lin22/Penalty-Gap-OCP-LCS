function FuncObj = create_FuncObj(self)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
% cost function
FuncObj.L_S = Function('L_S', {self.x, self.u, self.lambda}, {self.L_S}, {'x', 'u', 'lambda'}, {'L_S'});
% linear complementarity system
FuncObj.f = Function('f', {self.x, self.u, self.lambda}, {self.f}, {'x', 'u', 'lambda'}, {'f'});
FuncObj.g = Function('g', {self.x, self.u, self.lambda}, {self.g}, {'x', 'u', 'lambda'}, {'g'});
end