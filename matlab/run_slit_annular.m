function run_slit_annular()
% RUN_SLIT_ANNULAR Section 3.2 Slit Annular Plate Problem

addpath(genpath('.'));

fprintf('=== Section 3.2: Slit Annular Plate Problem ===\n');

Ri = 6.0; Ro = 10.0; a = 0.03;
E = 2.1e7; nu = 0.0;
P_max = 0.8;

nr = 5; ntheta = 40;
[X0, elements, Vn0] = mesh_slit_annular(Ri, Ro, nr, ntheta);

fixed_nodes = find(abs(atan2(X0(:, 2), X0(:, 1))) < 1e-4);
fixed_dofs = [];
for i = 1:length(fixed_nodes)
    fixed_dofs = [fixed_dofs, (fixed_nodes(i)-1)*5 + (1:5)];
end

theta_vals = atan2(X0(:, 2), X0(:, 1));
theta_vals(theta_vals < 0) = theta_vals(theta_vals < 0) + 2*pi;
loaded_nodes = find(abs(theta_vals - max(theta_vals)) < 1e-3);

load_dofs = [];
total_load = [];
for i = 1:length(loaded_nodes)
    load_dofs = [load_dofs; (loaded_nodes(i)-1)*5 + 3];
    total_load = [total_load; P_max / length(loaded_nodes)];
end

num_steps = 10;

fprintf('Running MITC4 (Slit Annular)...\n');
hist_m4 = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, false, Vn0);

fprintf('Running MITC4+ (Slit Annular)...\n');
hist_m4p = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, true, Vn0);

node_A = loaded_nodes(end);
hist_m4.disp  = abs(hist_m4.U((node_A-1)*5 + 3, :));
hist_m4p.disp = abs(hist_m4p.U((node_A-1)*5 + 3, :));

plot_results(hist_m4, hist_m4p, [], 'Slit Annular Plate Load-Displacement', 'slit_annular');

fprintf('Section 3.2 completed successfully.\n\n');

end
