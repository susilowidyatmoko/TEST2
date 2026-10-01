function run_scordelis()
% RUN_SCORDELIS Section 3.4 Scordelis-Lo Roof Problem

addpath(genpath('.'));

fprintf('=== Section 3.4: Scordelis-Lo Roof Problem ===\n');

R = 25.0; L = 50.0; theta_deg = 40; a = 0.25;
E = 4.32e8; nu = 0.0;
q_max = 4.0e4;

nx = 16; ny = 16;
[X0, elements, Vn0] = mesh_scordelis(R, L, theta_deg, nx, ny);

end_nodes = find(abs(X0(:, 1) - L/2) < 1e-5);
sym_x0 = find(abs(X0(:, 1)) < 1e-5);
sym_y0 = find(abs(X0(:, 2)) < 1e-5);

fixed_dofs = [];
for i = 1:length(end_nodes)
    fixed_dofs = [fixed_dofs, (end_nodes(i)-1)*5 + (1:3)];
end
for i = 1:length(sym_x0)
    fixed_dofs = [fixed_dofs, (sym_x0(i)-1)*5 + 1];
end
for i = 1:length(sym_y0)
    fixed_dofs = [fixed_dofs, (sym_y0(i)-1)*5 + 2];
end
fixed_dofs = unique(fixed_dofs);

Nnodes = size(X0, 1);
total_area = (L/2) * (R * theta_deg * pi / 180);
nodal_force = -q_max * total_area / Nnodes;

load_dofs = (1:Nnodes)*5 - 2;
total_load = nodal_force * ones(Nnodes, 1);

num_steps = 10;

fprintf('Running MITC4 (Scordelis-Lo Roof)...\n');
hist_m4 = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, false, Vn0);

fprintf('Running MITC4+ (Scordelis-Lo Roof)...\n');
hist_m4p = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, true, Vn0);

free_edge_nodes = find(abs(X0(:, 2) - max(X0(:, 2))) < 1e-4);
[~, idx_C] = min(abs(X0(free_edge_nodes, 1)));
node_C = free_edge_nodes(idx_C);

hist_m4.disp  = abs(hist_m4.U((node_C-1)*5 + 3, :));
hist_m4p.disp = abs(hist_m4p.U((node_C-1)*5 + 3, :));

plot_results(hist_m4, hist_m4p, [], 'Scordelis-Lo Roof Load-Displacement', 'scordelis');

fprintf('Section 3.4 completed successfully.\n\n');

end
