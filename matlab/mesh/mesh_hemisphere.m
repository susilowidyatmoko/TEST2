function [X, elements, Vn0] = mesh_hemisphere(R, phi_hole_deg, n_phi, n_theta)
% MESH_HEMISPHERE Generates mesh and exact spherical unit normal directors
% Inputs:
%   R            - Radius of hemisphere
%   phi_hole_deg - Angle of top hole in degrees (e.g. 18 deg)
%   n_phi, n_theta- Number of elements in phi and theta directions
% Outputs:
%   X            - Nnodes x 3 matrix
%   elements     - Nelems x 4 matrix
%   Vn0          - Nnodes x 3 unit normal directors

phi_min = phi_hole_deg * pi / 180;
phi_max = pi / 2;
theta_min = 0;
theta_max = pi / 2;

phi = linspace(phi_min, phi_max, n_phi + 1);
theta = linspace(theta_min, theta_max, n_theta + 1);

[TT, PP] = meshgrid(theta, phi);

XX = R * sin(PP) .* cos(TT);
YY = R * sin(PP) .* sin(TT);
ZZ = R * cos(PP);

X = [XX(:), YY(:), ZZ(:)];

elements = zeros(n_phi * n_theta, 4);
elem_idx = 1;

for j = 1:n_theta
    for i = 1:n_phi
        n1 = (j - 1) * (n_phi + 1) + i;
        n2 = (j - 1) * (n_phi + 1) + i + 1;
        n3 = j * (n_phi + 1) + i + 1;
        n4 = j * (n_phi + 1) + i;
        elements(elem_idx, :) = [n1, n2, n3, n4];
        elem_idx = elem_idx + 1;
    end
end

% Unit outward radial normal director for sphere
Vn0 = zeros(size(X));
for i = 1:size(X, 1)
    norm_x = norm(X(i, :));
    if norm_x > 1e-12
        Vn0(i, :) = X(i, :) / norm_x;
    else
        Vn0(i, :) = [0, 0, 1];
    end
end

end
