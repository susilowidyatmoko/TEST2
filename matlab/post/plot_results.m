function plot_results(history_mitc4, history_mitc4_plus, ref_data, title_str, filename_prefix)
% PLOT_RESULTS Plots load-displacement curves and deformed meshes for comparison
%
% Inputs:
%   history_mitc4      - Output struct from nonlinear_solver for MITC4
%   history_mitc4_plus - Output struct from nonlinear_solver for MITC4+
%   ref_data           - [disp, load] matrix for reference solution (optional)
%   title_str          - Title string for plots
%   filename_prefix    - File name prefix to save figures

if nargin < 3, ref_data = []; end
if nargin < 4, title_str = 'Load-Displacement Response'; end
if nargin < 5, filename_prefix = 'results'; end

% Figure 1: Load-Displacement Curve
fig1 = figure('Visible', 'off');
hold on; grid on; box on;

if ~isempty(ref_data)
    plot(ref_data(:, 1), ref_data(:, 2), 'k--', 'LineWidth', 2, 'DisplayName', 'Reference');
end

if ~isempty(history_mitc4)
    disp_m4 = history_mitc4.disp;
    load_m4 = history_mitc4.loads;
    plot(disp_m4, load_m4, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'MITC4');
end

if ~isempty(history_mitc4_plus)
    disp_m4p = history_mitc4_plus.disp;
    load_m4p = history_mitc4_plus.loads;
    plot(disp_m4p, load_m4p, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'MITC4+');
end

xlabel('Displacement');
ylabel('Load');
title(title_str);
legend('Location', 'best');

saveas(fig1, [filename_prefix, '_load_disp.png']);
close(fig1);

fprintf('Saved plot to %s_load_disp.png\n', filename_prefix);

end
