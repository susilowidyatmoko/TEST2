function run_twisted_cantilever()
% RUN_TWISTED_CANTILEVER Section 3.6 Twisted Cantilever Beam Problem

addpath(genpath('.'));

fprintf('=== Section 3.6: Twisted Cantilever Beam Problem ===\n');

L = 12.0; w = 1.0; a = 0.05;
total_twist_deg = 90;
E = 29.0e6; nu = 0.22;
P_max = 1.0;

nx = 12; ny = 2;
[X0, elements, Vn0] = mesh_twisted_cantilever(L, w, total_twist_deg, nx, ny);

root_nodes = find(abs(X0(:, 1)) < 1e-5);
fixed_dofs = [];
for i = 1:length(root_nodes)
    fixed_dofs = [fixed_dofs, (root_nodes(i)-1)*5 + (1:5)];
end

tip_nodes = find(abs(X0(:, 1) - L) < 1e-5);
load_dofs = [];
total_load = [];
for i = 1:length(tip_nodes)
    load_dofs = [load_dofs; (tip_nodes(i)-1)*5 + 3];
    total_load = [total_load; P_max / length(tip_nodes)];
end

num_steps = 10;

fprintf('Running MITC4 (Twisted Cantilever)...\n');
hist_m4 = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, false, Vn0);

fprintf('Running MITC4+ (Twisted Cantilever)...\n');
hist_m4p = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, true, Vn0);

tip_node_1 = tip_nodes(1);
hist_m4.disp  = abs(hist_m4.U((tip_node_1-1)*5 + 3, :));
hist_m4p.disp = abs(hist_m4p.U((tip_node_1-1)*5 + 3, :));

plot_results(hist_m4, hist_m4p, [], 'Twisted Cantilever Beam Load-Displacement', 'twisted_cantilever');

fprintf('Section 3.6 completed successfully.\n\n');

end
