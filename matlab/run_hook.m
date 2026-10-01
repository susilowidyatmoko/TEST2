function run_hook()
% RUN_HOOK Section 3.7 Hook Problem (Raasch Challenge)

addpath(genpath('.'));

fprintf('=== Section 3.7: Hook Problem (Raasch Challenge) ===\n');

width = 1.0; a = 0.05;
E = 3.3e3; nu = 0.35;
P_max = 100.0;

s_total = 14.12 * (150*pi/180) + 4.6 * (180*pi/180);
n_s = 20; n_w = 4;

[X0, elements, Vn0] = mesh_hook(width, s_total, n_s, n_w);

root_nodes = find(abs(X0(:, 1)) < 1e-4 & abs(X0(:, 2)) < 1e-4);
if isempty(root_nodes)
    root_nodes = 1:(n_w+1);
end

fixed_dofs = [];
for i = 1:length(root_nodes)
    fixed_dofs = [fixed_dofs, (root_nodes(i)-1)*5 + (1:5)];
end

tip_nodes = find(abs(X0(:, 1) - max(X0(:, 1))) < 1e-2);
load_dofs = [];
total_load = [];
for i = 1:length(tip_nodes)
    load_dofs = [load_dofs; (tip_nodes(i)-1)*5 + 3];
    total_load = [total_load; P_max / length(tip_nodes)];
end

num_steps = 10;

fprintf('Running MITC4 (Hook Problem)...\n');
hist_m4 = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, false, Vn0);

fprintf('Running MITC4+ (Hook Problem)...\n');
hist_m4p = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, true, Vn0);

tip_node_1 = tip_nodes(1);
hist_m4.disp  = abs(hist_m4.U((tip_node_1-1)*5 + 3, :));
hist_m4p.disp = abs(hist_m4p.U((tip_node_1-1)*5 + 3, :));

plot_results(hist_m4, hist_m4p, [], 'Hook Problem (Raasch) Load-Displacement', 'hook');

fprintf('Section 3.7 completed successfully.\n\n');

end
