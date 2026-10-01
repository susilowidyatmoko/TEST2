function run_cylinder()
% RUN_CYLINDER Section 3.5 Bending of Cylindrical Shell Structure

addpath(genpath('.'));

fprintf('=== Section 3.5: Cylindrical Shell Bending Problem ===\n');

R = 1.0; L = 1.0; theta_deg = 90;
E = 1.0e7; nu = 0.3;
M_max = 4.0e4;

aspect_ratios = [1/100, 1/1000];

for k = 1:length(aspect_ratios)
    a = R * aspect_ratios(k);
    fprintf('Testing thickness a/R = 1/%d ...\n', round(1/aspect_ratios(k)));

    [X0, elements, Vn0] = mesh_cylinder(R, L, theta_deg, 12, 12, false);
    fixed_nodes = find(abs(X0(:, 1)) < 1e-5);
    fixed_dofs = [];
    for i = 1:length(fixed_nodes)
        fixed_dofs = [fixed_dofs, (fixed_nodes(i)-1)*5 + (1:5)];
    end
    top_nodes = find(abs(X0(:, 1) - L) < 1e-5);
    load_dofs = [];
    total_load = [];
    for i = 1:length(top_nodes)
        load_dofs = [load_dofs; (top_nodes(i)-1)*5 + 4];
        total_load = [total_load; M_max / length(top_nodes)];
    end

    num_steps = 10;

    fprintf('  Running MITC4 (Uniform Mesh)...\n');
    hist_m4_uni = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, false, Vn0);

    fprintf('  Running MITC4+ (Uniform Mesh)...\n');
    hist_m4p_uni = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps, true, Vn0);

    [X0_dist, elements_dist, Vn0_dist] = mesh_cylinder(R, L, theta_deg, 12, 12, true);
    fixed_nodes_dist = find(abs(X0_dist(:, 1)) < 1e-5);
    fixed_dofs_dist = [];
    for i = 1:length(fixed_nodes_dist)
        fixed_dofs_dist = [fixed_dofs_dist, (fixed_nodes_dist(i)-1)*5 + (1:5)];
    end
    top_nodes_dist = find(abs(X0_dist(:, 1) - L) < 1e-5);
    load_dofs_dist = [];
    total_load_dist = [];
    for i = 1:length(top_nodes_dist)
        load_dofs_dist = [load_dofs_dist; (top_nodes_dist(i)-1)*5 + 4];
        total_load_dist = [total_load_dist; M_max / length(top_nodes_dist)];
    end

    fprintf('  Running MITC4 (Distorted Mesh)...\n');
    hist_m4_dist = nonlinear_solver(X0_dist, elements_dist, a, E, nu, fixed_dofs_dist, load_dofs_dist, total_load_dist, num_steps, false, Vn0_dist);

    fprintf('  Running MITC4+ (Distorted Mesh)...\n');
    hist_m4p_dist = nonlinear_solver(X0_dist, elements_dist, a, E, nu, fixed_dofs_dist, load_dofs_dist, total_load_dist, num_steps, true, Vn0_dist);

    top_node_1 = top_nodes(1);
    hist_m4_uni.disp  = abs(hist_m4_uni.U((top_node_1-1)*5 + 3, :));
    hist_m4p_uni.disp = abs(hist_m4p_uni.U((top_node_1-1)*5 + 3, :));

    plot_results(hist_m4_uni, hist_m4p_uni, [], sprintf('Cylinder Bending (a/R=1/%d)', round(1/aspect_ratios(k))), sprintf('cylinder_aR_%d', round(1/aspect_ratios(k))));
end

fprintf('Section 3.5 completed successfully.\n\n');

end
