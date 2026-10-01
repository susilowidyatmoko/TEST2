function [X, elements, Vn0] = mesh_hook(width, s_total, n_s, n_w)
% MESH_HOOK Generates mesh and unit normal directors for Hook problem (Raasch challenge)
% Inputs:
%   width   - Width along Z direction
%   s_total - Arc coordinate path
%   n_s, n_w- Number of elements along length and width
% Outputs:
%   X       - Nnodes x 3 matrix
%   elements- Nelems x 4 matrix
%   Vn0     - Nnodes x 3 unit normal directors

R1 = 14.12; phi1 = 150 * pi / 180;
R2 = 4.6;   phi2 = 180 * pi / 180;

L1 = R1 * phi1;
L2 = R2 * phi2;
L_tot = L1 + L2;

s = linspace(0, L_tot, n_s + 1);
z = linspace(-width/2, width/2, n_w + 1);

X_pts = zeros(n_s + 1, 3);
Vn_pts = zeros(n_s + 1, 3);

for i = 1:n_s + 1
    si = s(i);
    if si <= L1
        angle = si / R1;
        X_pts(i, 1) = R1 * (1 - cos(angle));
        X_pts(i, 2) = R1 * sin(angle);
        Vn_pts(i, :) = [-cos(angle), sin(angle), 0];
    else
        s2 = si - L1;
        angle2 = s2 / R2;
        C1 = [R1 * (1 - cos(phi1)), R1 * sin(phi1)];
        X_pts(i, 1) = C1(1) + R2 * (sin(phi1 - angle2) - sin(phi1));
        X_pts(i, 2) = C1(2) - R2 * (cos(phi1 - angle2) - cos(phi1));
        Vn_pts(i, :) = [cos(phi1 - angle2), sin(phi1 - angle2), 0];
    end
end

[ZZ, XX_line] = meshgrid(z, X_pts(:, 1));
[~,  YY_line] = meshgrid(z, X_pts(:, 2));
[~,  VnX_line]= meshgrid(z, Vn_pts(:, 1));
[~,  VnY_line]= meshgrid(z, Vn_pts(:, 2));
[~,  VnZ_line]= meshgrid(z, Vn_pts(:, 3));

X = [XX_line(:), YY_line(:), ZZ(:)];
Vn0 = [VnX_line(:), VnY_line(:), VnZ_line(:)];

elements = zeros(n_s * n_w, 4);
elem_idx = 1;

for j = 1:n_w
    for i = 1:n_s
        n1 = (j - 1) * (n_s + 1) + i;
        n2 = i + (j - 1) * (n_s + 1) + 1;
        n3 = i + j * (n_s + 1) + 1;
        n4 = i + j * (n_s + 1);
        elements(elem_idx, :) = [n1, n2, n3, n4];
        elem_idx = elem_idx + 1;
    end
end

end
