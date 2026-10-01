function run_cantilever()
% RUN_CANTILEVER Section 3.1 Cantilever Problem
% Evaluates both tip shearing force (Case a) and tip moment (Case b)

addpath(genpath('.'));

fprintf('=== Section 3.1: Cantilever Problem ===\n');

% Case (a): Tip shearing force
L = 10.0; b = 1.0; a = 0.1;
E = 1.2e6; nu = 0.0;
P_max = 25.0;
nx = 16; ny = 1;

[X0, elements, Vn0] = mesh_cantilever(L, b, nx, ny);

root_nodes = find(abs(X0(:, 1)) < 1e-6);
fixed_dofs = [];
for i = 1:length(root_nodes)
    fixed_dofs = [fixed_dofs, (root_nodes(i)-1)*5 + (1:5)];
end

tip_nodes = find(abs(X0(:, 1) - L) < 1e-6);
load_dofs_shear = [];
total_load_shear = [];
for i = 1:length(tip_nodes)
    load_dofs_shear = [load_dofs_shear; (tip_nodes(i)-1)*5 + 3]; % Z translation
    total_load_shear = [total_load_shear; P_max / length(tip_nodes)];
end

num_steps = 10;

fprintf('Running MITC4 (Tip Shear Force)...\n');
hist_m4_shear = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs_shear, total_load_shear, num_steps, false, Vn0);

fprintf('Running MITC4+ (Tip Shear Force)...\n');
hist_m4p_shear = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs_shear, total_load_shear, num_steps, true, Vn0);

tip_node_1 = tip_nodes(1);
hist_m4_shear.disp  = abs(hist_m4_shear.U((tip_node_1-1)*5 + 3, :));
hist_m4p_shear.disp = abs(hist_m4p_shear.U((tip_node_1-1)*5 + 3, :));

plot_results(hist_m4_shear, hist_m4p_shear, [], 'Cantilever Tip Shearing Force', 'cantilever_shear');

% Case (b): Tip moment
L_m = 12.0;
[X0_m, elements_m, Vn0_m] = mesh_cantilever(L_m, b, nx, ny);

root_nodes_m = find(abs(X0_m(:, 1)) < 1e-6);
fixed_dofs_m = [];
for i = 1:length(root_nodes_m)
    fixed_dofs_m = [fixed_dofs_m, (root_nodes_m(i)-1)*5 + (1:5)];
end

tip_nodes_m = find(abs(X0_m(:, 1) - L_m) < 1e-6);
load_dofs_m = [];
total_load_m = [];
I_val = b * a^3 / 12;
M_max = 2 * pi * E * I_val / L_m;

for i = 1:length(tip_nodes_m)
    load_dofs_m = [load_dofs_m; (tip_nodes_m(i)-1)*5 + 5]; % Beta rotation (Y axis)
    total_load_m = [total_load_m; M_max / length(tip_nodes_m)];
end

fprintf('Running MITC4 (Tip Moment)...\n');
hist_m4_mom = nonlinear_solver(X0_m, elements_m, a, E, nu, fixed_dofs_m, load_dofs_m, total_load_m, num_steps, false, Vn0_m);

fprintf('Running MITC4+ (Tip Moment)...\n');
hist_m4p_mom = nonlinear_solver(X0_m, elements_m, a, E, nu, fixed_dofs_m, load_dofs_m, total_load_m, num_steps, true, Vn0_m);

tip_node_m1 = tip_nodes_m(1);
hist_m4_mom.disp  = abs(hist_m4_mom.U((tip_node_m1-1)*5 + 3, :));
hist_m4p_mom.disp = abs(hist_m4p_mom.U((tip_node_m1-1)*5 + 3, :));

plot_results(hist_m4_mom, hist_m4p_mom, [], 'Cantilever Tip Moment', 'cantilever_moment');

fprintf('Section 3.1 completed successfully.\n\n');

end
