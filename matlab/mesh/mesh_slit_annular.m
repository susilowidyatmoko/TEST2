function [X, elements, Vn0] = mesh_slit_annular(Ri, Ro, nr, ntheta)
% MESH_SLIT_ANNULAR Generates mesh and unit normal directors for slit annular plate
% Inputs:
%   Ri, Ro       - Inner and outer radii
%   nr, ntheta   - Number of elements along radial and circumferential directions
% Outputs:
%   X            - Nnodes x 3 matrix
%   elements     - Nelems x 4 matrix
%   Vn0          - Nnodes x 3 unit normal directors

r = linspace(Ri, Ro, nr + 1);
theta = linspace(0, 2*pi, ntheta + 1);

[TT, RR] = meshgrid(theta, r);
XX = RR .* cos(TT);
YY = RR .* sin(TT);
ZZ = zeros(size(XX));

X = [XX(:), YY(:), ZZ(:)];

elements = zeros(nr * ntheta, 4);
elem_idx = 1;

for j = 1:ntheta
    for i = 1:nr
        n1 = (j - 1) * (nr + 1) + i;
        n2 = (j - 1) * (nr + 1) + i + 1;
        n3 = j * (nr + 1) + i + 1;
        n4 = j * (nr + 1) + i;
        elements(elem_idx, :) = [n1, n2, n3, n4];
        elem_idx = elem_idx + 1;
    end
end

Nnodes = size(X, 1);
Vn0 = repmat([0, 0, 1], Nnodes, 1);

end
