function [X_new, V_n_new, V_1_new, V_2_new] = update_geometry(X, V_n, V_1, V_2, U)
% UPDATE_GEOMETRY Update nodal coordinates and orthogonal director vectors
% Inputs:
%   X     - Nnodes x 3 matrix of nodal positions
%   V_n   - Nnodes x 3 unit normal directors
%   V_1   - Nnodes x 3 unit tangent directors 1
%   V_2   - Nnodes x 3 unit tangent directors 2
%   U     - 5*Nnodes x 1 displacement increment vector [u, v, w, alpha, beta]
% Outputs:
%   X_new   - Updated positions
%   V_n_new - Updated normal vectors
%   V_1_new - Updated tangent vectors 1
%   V_2_new - Updated tangent vectors 2

Nnodes = size(X, 1);
X_new = X;
V_n_new = V_n;
V_1_new = V_1;
V_2_new = V_2;

for i = 1:Nnodes
    idx = (i-1)*5 + (1:5);
    u_i = U(idx(1:3));
    alpha_i = U(idx(4));
    beta_i  = U(idx(5));

    % Update position
    X_new(i, :) = X(i, :) + u_i';

    % Rotation vector theta = -alpha * V2 + beta * V1
    theta_vec = -alpha_i * V_2(i, :)' + beta_i * V_1(i, :)';
    theta_mag = norm(theta_vec);

    if theta_mag > 1e-12
        axis_rot = theta_vec / theta_mag;
        q0 = cos(theta_mag / 2);
        q1 = axis_rot(1) * sin(theta_mag / 2);
        q2 = axis_rot(2) * sin(theta_mag / 2);
        q3 = axis_rot(3) * sin(theta_mag / 2);

        % Quaternion rotation matrix Q
        Q = 2 * [q0^2 + q1^2 - 0.5,  q1*q2 - q0*q3,    q3*q1 + q0*q2;
                 q1*q2 + q0*q3,    q0^2 + q2^2 - 0.5,  q2*q3 - q0*q1;
                 q3*q1 - q0*q2,    q2*q3 + q0*q1,    q0^2 + q3^2 - 0.5];

        % Rotate director vectors
        V_n_new(i, :) = (Q * V_n(i, :)')';
        V_1_new(i, :) = (Q * V_1(i, :)')';
        V_2_new(i, :) = (Q * V_2(i, :)')';
    end

    % Re-normalize director vectors to prevent numerical drift
    V_n_new(i, :) = V_n_new(i, :) / norm(V_n_new(i, :));

    % Orthogonalize V_1 relative to V_n
    V_1_proj = V_1_new(i, :)' - dot(V_1_new(i, :), V_n_new(i, :)) * V_n_new(i, :)';
    if norm(V_1_proj) > 1e-8
        V_1_new(i, :) = (V_1_proj / norm(V_1_proj))';
    else
        % Fallback if collinear
        if abs(V_n_new(i, 1)) < 0.8
            ref_v = [1, 0, 0];
        else
            ref_v = [0, 1, 0];
        end
        v1_tmp = cross(ref_v, V_n_new(i, :));
        V_1_new(i, :) = v1_tmp / norm(v1_tmp);
    end

    % Compute V_2 = V_n x V_1
    v2_tmp = cross(V_n_new(i, :), V_1_new(i, :));
    V_2_new(i, :) = v2_tmp / norm(v2_tmp);
end

end
