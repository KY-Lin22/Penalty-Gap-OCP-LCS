function [dz, gamma_h_k, Info] = evaluate_search_direction(self, h, J_grad, J_penalty_hessian)
%evaluate the search direction 

% KKT residual
KKT_residual = [-J_grad'; -h];
% KKT matrix
[i_J_pen_hess, j_J_pen_hess, s_J_pen_hess] = find(J_penalty_hessian);
KKT_matrix_update = sparse(i_J_pen_hess, j_J_pen_hess, s_J_pen_hess,...
    self.NLP.Dim.z + self.NLP.Dim.h, self.NLP.Dim.z + self.NLP.Dim.h, length(s_J_pen_hess));
KKT_matrix = self.KKT_matrix_constant + KKT_matrix_update;
% solve linear system
dz_gamma_h_k = KKT_matrix\KKT_residual;
dz = dz_gamma_h_k(1 : self.NLP.Dim.z, 1);
gamma_h_k = dz_gamma_h_k(self.NLP.Dim.z + 1 : end, 1);

% create info
Info.status = 1;
Info.msg = 'succeeds';

end