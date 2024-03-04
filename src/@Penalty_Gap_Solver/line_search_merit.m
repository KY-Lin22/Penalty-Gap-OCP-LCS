function [z_k, Info] = line_search_merit(self, beta, z, dz, p, J, h, J_grad, J_penalty_hessian)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% some quantities at current iterate z
% directional derivative of cost
J_DD = J_grad * dz;
% constraint violation M 
% - L1 norm follows IPOPT, and also the cost is the sum of stage cost
if self.Option.LineSearch.scaling_constraint_violation
    % - scaled by time step to consistent with the ocp cost that has been scaled
    M = self.OCP.timeStep * norm(h, 1);
else
    M = norm(h, 1);
end
% penalty parameter
if self.Option.LineSearch.using_Hessian
    J_ocp_hessian = self.NLP.FuncObj.J_ocp_hessian;
    beta_Trial = (J_DD + 1/2*dz'*(J_ocp_hessian + J_penalty_hessian)*dz)/((1 - self.Option.LineSearch.rho) * M);
else
    beta_Trial = (J_DD)/((1 - self.Option.LineSearch.rho) * M);
end

if beta >= beta_Trial
    beta_k = beta;
else
    beta_k = beta_Trial + 1;
end
% merit and its directional derivative 
merit = J + beta_k * M;
merit_DD = J_DD - beta_k * M;

%% backtracking line search
has_found_new_iterate = false;
stepSize_init = 1;

while ~has_found_new_iterate
     %% Step 1: estimate trial stepsize, iterate, cost, infeasibility and merit
     % step size and z
     stepSize_trial = max([stepSize_init, self.Option.LineSearch.stepSize_Min]);
     z_trial = z + stepSize_trial * dz;
     % cost and constraint
     J_trial = full(self.NLP.FuncObj.J(z_trial, p));
     h_trial = full(self.NLP.FuncObj.h(z_trial, p));
     % constraint infeasibility
     if self.Option.LineSearch.scaling_constraint_violation
         M_trial = self.OCP.timeStep * norm(h_trial, 1);
     else
         M_trial = norm(h_trial, 1);
     end
     % merit
     merit_trial = J_trial + beta_k * M_trial;

     %% Step 2: check sufficient decrease condition
     if merit_trial <= merit + stepSize_trial * self.Option.LineSearch.nu_D * merit_DD
         has_found_new_iterate = true;
         status = 1;
     end

     %% Step 3: checking min stepsize
    if ~has_found_new_iterate
        if stepSize_trial == self.Option.LineSearch.stepSize_Min
            % linesearch fails on the min stepsize, break backtracking linesearch procedure
            status = 0;
            break
        else
            % estimate a smaller stepsize
            stepSize_init = self.Option.LineSearch.stepSize_DecayRate * stepSize_init;
        end
    end     

end

%% organize output
Info.status = status;
switch status
    case 0
        % fail, return the previous one
        z_k = z;        
    case 1
        % success, return the new iterate
        z_k = z_trial;
        Info.J = J_trial;
        Info.h = h_trial;
        Info.beta = beta_k;
        Info.stepSize = stepSize_trial;
        Info.merit = [merit, merit_trial];
end

end