function [X, elements, Vn0] = mesh_scordelis(R, L, theta_deg, nx, ny)
% MESH_SCORDELIS Generates 1/4 mesh and exact unit normal directors for Scordelis-Lo roof
% Inputs:
%   R         - Cylinder radius
%   L         - Total length (1/2 modeled due to symmetry)
%   theta_deg - Semi-angle of roof in degrees (e.g. 40 deg)
%   nx, ny    - Elements along axial and circumferential directions
% Outputs:
%   X         - Nnodes x 3 matrix
%   elements  - Nelems x 4 matrix
%   Vn0       - Nnodes x 3 unit normal directors

theta_rad = theta_deg * pi / 180;

x = linspace(0, L/2, nx + 1);
theta = linspace(0, theta_rad, ny + 1);

[TT, XX] = meshgrid(theta, x);

YY = R * sin(TT);
ZZ = R * cos(TT);

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

% Radial unit normal director for cylinder
Vn0 = zeros(size(X));
TT_vec = TT(:);
for i = 1:size(X, 1)
    Vn0(i, :) = [0, sin(TT_vec(i)), cos(TT_vec(i))];
end

end
