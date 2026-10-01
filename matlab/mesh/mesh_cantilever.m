function [X, elements, Vn0] = mesh_cantilever(L, b, nx, ny)
% MESH_CANTILEVER Generates mesh and unit normal directors for cantilever beam
% Inputs:
%   L, b   - Length (X) and Width (Y)
%   nx, ny - Number of elements along length and width
% Outputs:
%   X        - (nx+1)*(ny+1) x 3 nodal coordinates matrix
%   elements - nx*ny x 4 connectivity matrix
%   Vn0      - Nnodes x 3 unit normal directors

x = linspace(0, L, nx + 1);
y = linspace(-b/2, b/2, ny + 1);

[YY, XX] = meshgrid(y, x);
X = [XX(:), YY(:), zeros(numel(XX), 1)];

elements = zeros(nx * ny, 4);
elem_idx = 1;

for i = 1:nx
    for j = 1:ny
        n1 = (i - 1) * (ny + 1) + j;
        n2 = i * (ny + 1) + j;
        n3 = i * (ny + 1) + j + 1;
        n4 = (i - 1) * (ny + 1) + j + 1;
        elements(elem_idx, :) = [n1, n2, n3, n4];
        elem_idx = elem_idx + 1;
    end
end

Nnodes = size(X, 1);
Vn0 = repmat([0, 0, 1], Nnodes, 1);

end
