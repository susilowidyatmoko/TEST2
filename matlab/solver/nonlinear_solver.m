function [history] = nonlinear_solver(X0, elements, thick, E, nu, fixed_dofs, load_dofs, total_load, num_steps, is_mitc4_plus, Vn0, max_iter, tol)
% NONLINEAR_SOLVER Geometrically nonlinear Total Lagrangian solver using Newton-Raphson
%
% Inputs:
%   X0           - Nnodes x 3 initial nodal coordinates
%   elements     - Nelems x 4 node connectivity matrix
%   thick        - Nnodes x 1 thickness array
%   E            - Young's modulus
%   nu           - Poisson's ratio
%   fixed_dofs   - Vector of global DOFs to fix (1-based, 5 DOFs per node)
%   load_dofs    - Vector of global DOFs where loads are applied
%   total_load   - Vector of total applied load at load_dofs
%   num_steps    - Number of load increments
%   is_mitc4_plus- boolean (true for MITC4+, false for MITC4)
%   Vn0          - Nnodes x 3 initial unit normal directors (optional)
%   max_iter     - Maximum Newton-Raphson iterations per step (default 25)
%   tol          - Convergence tolerance (default 1e-4)

if nargin < 10, is_mitc4_plus = true; end
if nargin < 11 || isempty(Vn0), Vn0 = repmat([0, 0, 1], size(X0, 1), 1); end
if nargin < 12, max_iter = 25; end
if nargin < 13, tol = 1e-4; end

Nnodes = size(X0, 1);
Nelems = size(elements, 1);
Ndofs = 5 * Nnodes;

if numel(thick) == 1
    thick = thick * ones(Nnodes, 1);
end

% Compute initial orthogonal tangent directors V10, V20
V10 = zeros(Nnodes, 3);
V20 = zeros(Nnodes, 3);

for i = 1:Nnodes
    Vn0(i, :) = Vn0(i, :) / norm(Vn0(i, :));
    if abs(Vn0(i, 1)) < 0.8
        ref_v = [1, 0, 0];
    else
        ref_v = [0, 1, 0];
    end
    v1_tmp = cross(ref_v, Vn0(i, :));
    V10(i, :) = v1_tmp / norm(v1_tmp);
    v2_tmp = cross(Vn0(i, :), V10(i, :));
    V20(i, :) = v2_tmp / norm(v2_tmp);
end

C_mat = material_matrix(E, nu);

Xt = X0;
Vnt = Vn0;
V1t = V10;
V2t = V20;
U_total = zeros(Ndofs, 1);

all_dofs = 1:Ndofs;
free_dofs = setdiff(all_dofs, fixed_dofs);

history.steps = zeros(num_steps, 1);
history.U = zeros(Ndofs, num_steps);
history.Xt = cell(num_steps, 1);
history.loads = zeros(length(load_dofs), num_steps);

for step = 1:num_steps
    step_factor = step / num_steps;
    F_ext = zeros(Ndofs, 1);
    F_ext(load_dofs) = step_factor * total_load;

    iter = 0;
    converged = false;
    rel_res = 1.0;

    while ~converged && iter < max_iter
        iter = iter + 1;

        K_glob = zeros(Ndofs, Ndofs);
        F_int = zeros(Ndofs, 1);

        for e = 1:Nelems
            nodes_e = elements(e, :);
            dof_e = zeros(1, 20);
            for n = 1:4
                dof_e((n-1)*5 + (1:5)) = (nodes_e(n)-1)*5 + (1:5);
            end

            X0_e  = X0(nodes_e, :);
            Xt_e  = Xt(nodes_e, :);
            Vn0_e = Vn0(nodes_e, :); V10_e = V10(nodes_e, :); V20_e = V20(nodes_e, :);
            Vnt_e = Vnt(nodes_e, :); V1t_e = V1t(nodes_e, :); V2t_e = V2t(nodes_e, :);
            thick_e = thick(nodes_e);

            [Ke, Fe] = mitc4_plus_element(X0_e, Xt_e, Vn0_e, V10_e, V20_e, Vnt_e, V1t_e, V2t_e, thick_e, C_mat, is_mitc4_plus);

            K_glob(dof_e, dof_e) = K_glob(dof_e, dof_e) + Ke;
            F_int(dof_e) = F_int(dof_e) + Fe;
        end

        Residual = F_ext - F_int;
        res_norm = norm(Residual(free_dofs));
        ext_norm = norm(F_ext(free_dofs));
        if ext_norm < 1e-10, ext_norm = 1.0; end

        rel_res = res_norm / ext_norm;

        if rel_res < tol && iter > 1
            converged = true;
            break;
        end

        dU_free = K_glob(free_dofs, free_dofs) \ Residual(free_dofs);
        dU = zeros(Ndofs, 1);
        dU(free_dofs) = dU_free;

        [Xt, Vnt, V1t, V2t] = update_geometry(Xt, Vnt, V1t, V2t, dU);
        U_total = U_total + dU;
    end

    fprintf('Step %2d/%2d: Iterations = %2d, Relative Residual = %.4e\n', step, num_steps, iter, rel_res);

    history.steps(step) = step;
    history.U(:, step) = U_total;
    history.Xt{step} = Xt;
    history.loads(:, step) = F_ext(load_dofs);
end

end
