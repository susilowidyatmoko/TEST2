function C = material_matrix(E, nu, k_shear)
% MATERIAL_MATRIX Constitutive matrix for 3D shell in local coordinates
% Inputs:
%   E       - Young's modulus
%   nu      - Poisson's ratio
%   k_shear - Transverse shear factor (default: 5/6)
% Output:
%   C       - 5x5 material constitutive matrix

if nargin < 3
    k_shear = 5/6;
end

% Local strain vector order: [e11, e22, 2*e12, 2*e13, 2*e23]^T
C = zeros(5, 5);
factor = E / (1 - nu^2);

% In-plane membrane and bending terms
C(1,1) = factor;
C(1,2) = factor * nu;
C(2,1) = factor * nu;
C(2,2) = factor;
C(3,3) = factor * (1 - nu) / 2;

% Transverse shear terms
G = E / (2 * (1 + nu));
C(4,4) = k_shear * G;
C(5,5) = k_shear * G;

end
