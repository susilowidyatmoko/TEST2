function [X, elements, Vn0] = mesh_cylinder(R, L, theta_deg, nx, ny, is_distorted)
% MESH_CYLINDER Generates cylinder mesh and unit normal directors
% Inputs:
%   R, L, theta_deg - Radius, length, arc angle
%   nx, ny          - Elements along x and y directions
%   is_distorted    - Boolean flag for distorted mesh pattern
% Outputs:
%   X               - Nnodes x 3 matrix
%   elements        - Nelems x 4 matrix
%   Vn0             - Nnodes x 3 unit normal directors

if nargin < 6, is_distorted = false; end

theta_rad = theta_deg * pi / 180;

if ~is_distorted
    x = linspace(0, L, nx + 1);
    s = linspace(0, theta_rad, ny + 1);
else
    ratio_x = 1:nx;
    cumsum_x = [0, cumsum(ratio_x)] / sum(ratio_x);
    x = L * cumsum_x;

    ratio_s = 1:ny;
    cumsum_s = [0, cumsum(ratio_s)] / sum(ratio_s);
    s = theta_rad * cumsum_s;
end

[SS, XX] = meshgrid(s, x);

YY = R * sin(SS);
ZZ = R * cos(SS);

X = [XX(:), YY(:), ZZ(:)];

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
SS_vec = SS(:);
for i = 1:size(X, 1)
    Vn0(i, :) = [0, sin(SS_vec(i)), cos(SS_vec(i))];
end

end
