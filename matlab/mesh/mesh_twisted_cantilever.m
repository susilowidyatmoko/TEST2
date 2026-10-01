function [X, elements, Vn0] = mesh_twisted_cantilever(L, w, total_twist_deg, nx, ny)
% MESH_TWISTED_CANTILEVER Generates twisted cantilever mesh and exact unit normal directors
% Inputs:
%   L, w            - Length along X, width along Y before twist
%   total_twist_deg - Total twist angle in degrees
%   nx, ny          - Number of elements along length and width
% Outputs:
%   X               - Nnodes x 3 matrix
%   elements        - Nelems x 4 matrix
%   Vn0             - Nnodes x 3 unit normal directors

x = linspace(0, L, nx + 1);
y = linspace(-w/2, w/2, ny + 1);

[YY, XX] = meshgrid(y, x);

theta_x = (total_twist_deg * pi / 180) * (XX / L);

XX_twist = XX;
YY_twist = YY .* cos(theta_x);
ZZ_twist = YY .* sin(theta_x);

X = [XX_twist(:), YY_twist(:), ZZ_twist(:)];

elements = zeros(nx * ny, 4);
elem_idx = 1;

for j = 1:ny
    for i = 1:nx
        n1 = (j - 1) * (nx + 1) + i;
        n2 = i + (j - 1) * (nx + 1) + 1;
        n3 = i + j * (nx + 1) + 1;
        n4 = i + j * (nx + 1);
        elements(elem_idx, :) = [n1, n2, n3, n4];
        elem_idx = elem_idx + 1;
    end
end

Vn0 = zeros(size(X));
theta_vec = theta_x(:);
for i = 1:size(X, 1)
    ang = theta_vec(i);
    Vn0(i, :) = [0, -sin(ang), cos(ang)];
end

end
