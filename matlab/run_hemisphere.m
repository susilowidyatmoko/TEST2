function run_hemisphere()
% RUN_HEMISPHERE Section 3.3 Hemispherical Shell with Open Hole

addpath(genpath('.'));

fprintf('=== Section 3.3: Hemispherical Shell Problem ===\n');

R = 10.0; a = 0.04;
E = 6.825e7; nu = 0.3;
P_max = 400.0;

n_phi = 16; n_theta = 16;
[X0, elements, Vn0] = mesh_hemisphere(R, 18, n_phi, n_theta);

edge_y0 = find(abs(X0(:, 2)) < 1e-5);
edge_x0 = find(abs(X0(:, 1)) < 1e-5);

fixed_dofs = [];
for i = 1:length(edge_y0)
    fixed_dofs = [fixed_dofs, (edge_y0(i)-1)*5 + 2];
end
for i = 1:length(edge_x0)
    fixed_dofs = [fixed_dofs, (edge_x0(i)-1)*5 + 1];
end
fixed_dofs = unique(fixed_dofs);

rim_nodes = find(abs(X0(:, 3)) < 1e-4);
[~, idx_A] = min(abs(X0(rim_nodes, 2)));
[~, idx_B] = min(abs(X0(rim_nodes, 1)));

node_A = rim_nodes(idx_A);
node_B = rim_nodes(idx_B);

load_dofs = [(node_A-1)*5 + 1; (node_B-1)*5 + 2];
total_load = [-P_max / 2; P_max / 2];

num_steps = 10;

fprintf('Running MITC4 (Hemisphere)...\n');
hist_m4 = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, false, Vn0);

fprintf('Running MITC4+ (Hemisphere)...\n');
hist_m4p = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, true, Vn0);

hist_m4.disp  = abs(hist_m4.U((node_A-1)*5 + 1, :));
hist_m4p.disp = abs(hist_m4p.U((node_A-1)*5 + 1, :));

plot_results(hist_m4, hist_m4p, [], 'Hemispherical Shell Load-Displacement', 'hemisphere');

fprintf('Section 3.3 completed successfully.\n\n');

end
