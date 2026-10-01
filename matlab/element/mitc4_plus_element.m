function [Ke, Fe, S_local] = mitc4_plus_element(X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, C_mat, is_mitc4_plus)
% MITC4_PLUS_ELEMENT Tangent stiffness matrix and internal force vector for
% MITC4 / MITC4+ geometrically nonlinear continuum mechanics-based shell element.
%
% Inputs:
%   X0           - 4x3 initial nodal coordinates
%   Xt           - 4x3 current nodal coordinates
%   Vn0, V10, V20- 4x3 initial director vectors at 4 nodes
%   Vnt, V1t, V2t- 4x3 current director vectors at 4 nodes
%   thick        - 4x1 or scalar nodal thicknesses
%   C_mat        - 5x5 material constitutive matrix
%   is_mitc4_plus- boolean flag: true for MITC4+, false for classic MITC4
%
% Outputs:
%   Ke           - 20x20 element tangent stiffness matrix
%   Fe           - 20x1 element internal force vector
%   S_local      - 2x2x2x5 local Piola-Kirchhoff stresses at Gauss points

if nargin < 11
    is_mitc4_plus = true;
end

if numel(thick) == 1
    thick = thick * ones(4, 1);
end

Ke = zeros(20, 20);
Fe = zeros(20, 1);

% 2x2x2 Gauss integration points and weights
gp = 1/sqrt(3) * [-1, 1];
gw = [1, 1];

% Characteristic vectors for MITC4+ membrane strain field evaluation
n_vec  = [-1;  1;  1; -1];
g_vec  = [-1; -1;  1;  1];
ng_vec = [ 1; -1;  1; -1];

xr_t = 0.25 * (Xt' * n_vec);   % 3x1
xs_t = 0.25 * (Xt' * g_vec);   % 3x1
xd_t = 0.25 * (Xt' * ng_vec);  % 3x1

xr_0 = 0.25 * (X0' * n_vec);
xs_0 = 0.25 * (X0' * g_vec);
xd_0 = 0.25 * (X0' * ng_vec);

% Plane P normal at configuration t
n_plane_t = cross(xr_t, xs_t);
norm_n_t = norm(n_plane_t);
if norm_n_t > 1e-12
    n_plane_t = n_plane_t / norm_n_t;
else
    n_plane_t = [0; 0; 1];
end

% Dual base vectors mr_t, ms_t on plane P
A_dual_t = [xr_t'; xs_t'; n_plane_t'];
mr_t = A_dual_t \ [1; 0; 0];
ms_t = A_dual_t \ [0; 1; 0];

cr_t = dot(mr_t, xd_t);
cs_t = dot(ms_t, xd_t);
cd_t = cr_t^2 + cs_t^2 - 1;

if abs(cd_t) < 1e-8
    cd_t = sign(cd_t + 1e-15) * 1e-8;
end

aA_t = cr_t * (cr_t - 1) / (2 * cd_t);
aB_t = cr_t * (cr_t + 1) / (2 * cd_t);
aC_t = cs_t * (cs_t - 1) / (2 * cd_t);
aD_t = cs_t * (cs_t + 1) / (2 * cd_t);
aE_t = 2 * cr_t * cs_t / cd_t;

% Same parameters for configuration at time 0
n_plane_0 = cross(xr_0, xs_0);
norm_n_0 = norm(n_plane_0);
if norm_n_0 > 1e-12
    n_plane_0 = n_plane_0 / norm_n_0;
else
    n_plane_0 = [0; 0; 1];
end
A_dual_0 = [xr_0'; xs_0'; n_plane_0'];
mr_0 = A_dual_0 \ [1; 0; 0];
ms_0 = A_dual_0 \ [0; 1; 0];

cr_0 = dot(mr_0, xd_0);
cs_0 = dot(ms_0, xd_0);
cd_0 = cr_0^2 + cs_0^2 - 1;
if abs(cd_0) < 1e-8
    cd_0 = sign(cd_0 + 1e-15) * 1e-8;
end

aA_0 = cr_0 * (cr_0 - 1) / (2 * cd_0);
aB_0 = cr_0 * (cr_0 + 1) / (2 * cd_0);
aC_0 = cs_0 * (cs_0 - 1) / (2 * cd_0);
aD_0 = cs_0 * (cs_0 + 1) / (2 * cd_0);
aE_0 = 2 * cr_0 * cs_0 / cd_0;

% Tying positions for transverse shear (A, B, C, D)
r_tying_shear = [0,  0, -1,  1];
s_tying_shear = [-1, 1,  0,  0];

% Tying positions for membrane strains (A, B, C, D, E)
r_tying_mb = [ 0, 0, -1, 1, 0];
s_tying_mb = [-1, 1,  0, 0, 0];

% Evaluate tying point quantities
B_shear_tied = cell(4, 1);
g_rt_tied_0 = zeros(4, 1); g_rt_tied_t = zeros(4, 1);
g_st_tied_0 = zeros(4, 1); g_st_tied_t = zeros(4, 1);

for i = 1:4
    r_i = r_tying_shear(i);
    s_i = s_tying_shear(i);
    [~, B_sh, g_cov_0, g_cov_t] = eval_kinematics(r_i, s_i, 0, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick);
    B_shear_tied{i} = B_sh;
    g_rt_tied_0(i) = g_cov_0(4);
    g_st_tied_0(i) = g_cov_0(5);
    g_rt_tied_t(i) = g_cov_t(4);
    g_st_tied_t(i) = g_cov_t(5);
end

B_mb_tied = cell(5, 1);
g_rr_mb_0 = zeros(5, 1); g_rr_mb_t = zeros(5, 1);
g_ss_mb_0 = zeros(5, 1); g_ss_mb_t = zeros(5, 1);
g_rs_mb_0 = zeros(5, 1); g_rs_mb_t = zeros(5, 1);

for i = 1:5
    r_i = r_tying_mb(i);
    s_i = s_tying_mb(i);
    [~, B_m, g_cov_0, g_cov_t] = eval_kinematics(r_i, s_i, 0, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick);
    B_mb_tied{i} = B_m(1:3, :);
    g_rr_mb_0(i) = g_cov_0(1);
    g_ss_mb_0(i) = g_cov_0(2);
    g_rs_mb_0(i) = g_cov_0(3);
    g_rr_mb_t(i) = g_cov_t(1);
    g_ss_mb_t(i) = g_cov_t(2);
    g_rs_mb_t(i) = g_cov_t(3);
end

S_local = zeros(2, 2, 2, 5);

% Integration over element volume (2x2x2 Gauss integration)
for ig = 1:2
    r = gp(ig);
    wr = gw(ig);
    for jg = 1:2
        s = gp(jg);
        ws = gw(jg);
        for kg = 1:2
            zeta = gp(kg);
            wz = gw(kg);

            [J0, B_direct, g_cov_0, g_cov_t, g_inv_0, L1, L2, L3] = eval_kinematics(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick);

            detJ0 = det(J0);
            dV = detJ0 * wr * ws * wz;

            B_cov = B_direct;

            % Transverse shear strain interpolation (MITC4)
            B_cov(4, :) = 0.5*(1+s)*B_shear_tied{1}(4, :) + 0.5*(1-s)*B_shear_tied{2}(4, :);
            B_cov(5, :) = 0.5*(1+r)*B_shear_tied{3}(5, :) + 0.5*(1-r)*B_shear_tied{4}(5, :);

            % Membrane strain interpolation (MITC4+)
            if is_mitc4_plus
                B_rr = 0.5*(1 - 2*aA_0 + s + 2*aA_0*s^2)*B_mb_tied{1}(1,:) ...
                     + 0.5*(1 - 2*aB_0 - s + 2*aB_0*s^2)*B_mb_tied{2}(1,:) ...
                     + aC_0*(-1 + s^2)*B_mb_tied{3}(2,:) ...
                     + aD_0*(-1 + s^2)*B_mb_tied{4}(2,:) ...
                     + aE_0*(-1 + s^2)*B_mb_tied{5}(3,:);

                B_ss = aA_0*(-1 + r^2)*B_mb_tied{1}(1,:) ...
                     + aB_0*(-1 + r^2)*B_mb_tied{2}(1,:) ...
                     + 0.5*(1 - 2*aC_0 + r + 2*aC_0*r^2)*B_mb_tied{3}(2,:) ...
                     + 0.5*(1 - 2*aD_0 - r + 2*aD_0*r^2)*B_mb_tied{4}(2,:) ...
                     + aE_0*(-1 + r^2)*B_mb_tied{5}(3,:);

                B_rs = 0.25*(r + 4*aA_0*r*s)*B_mb_tied{1}(1,:) ...
                     + 0.25*(-r + 4*aB_0*r*s)*B_mb_tied{2}(1,:) ...
                     + 0.25*(s + 4*aC_0*r*s)*B_mb_tied{3}(2,:) ...
                     + 0.25*(-s + 4*aD_0*r*s)*B_mb_tied{4}(2,:) ...
                     + (1 + aE_0*r*s)*B_mb_tied{5}(3,:);

                [~, B_m_direct] = eval_kinematics(r, s, 0, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick);
                B_cov(1:3, :) = B_cov(1:3, :) - B_m_direct(1:3, :) + [B_rr; B_ss; 2*B_rs];
            end

            T = compute_covariant_to_local_transform(L1, L2, L3, g_inv_0);
            B_local = T * B_cov;

            E_cov = compute_covariant_GL_strain(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, ...
                is_mitc4_plus, aA_0, aB_0, aC_0, aD_0, aE_0, aA_t, aB_t, aC_t, aD_t, aE_t, ...
                g_rt_tied_0, g_st_tied_0, g_rt_tied_t, g_st_tied_t, ...
                g_rr_mb_0, g_ss_mb_0, g_rs_mb_0, g_rr_mb_t, g_ss_mb_t, g_rs_mb_t);

            E_local = T * E_cov;

            S_vec = C_mat * E_local;
            S_local(ig, jg, kg, :) = S_vec;

            K_geo = compute_geometric_stiffness(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, S_vec, T, g_inv_0);

            Ke = Ke + (B_local' * C_mat * B_local + K_geo) * dV;
            Fe = Fe + (B_local' * S_vec) * dV;
        end
    end
end

end


function [J0, B, g_cov_0, g_cov_t, g_inv_0, L1, L2, L3] = eval_kinematics(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick)
% Evaluates shape functions, derivatives, Jacobians, covariant base vectors,
% local Cartesian axes, and strain-displacement matrix B.

h = 0.25 * [(1-r)*(1-s); (1+r)*(1-s); (1+r)*(1+s); (1-r)*(1+s)];
dh_dr = 0.25 * [-(1-s);  (1-s);  (1+s); -(1+s)];
dh_ds = 0.25 * [-(1-r); -(1+r);  (1+r);  (1-r)];

dx_m_0_dr = X0' * dh_dr;
dx_m_0_ds = X0' * dh_ds;

x_b_0 = 0.5 * (Vn0 .* thick)' * h;
dx_b_0_dr = 0.5 * (Vn0 .* thick)' * dh_dr;
dx_b_0_ds = 0.5 * (Vn0 .* thick)' * dh_ds;

g1_0 = dx_m_0_dr + zeta * dx_b_0_dr;
g2_0 = dx_m_0_ds + zeta * dx_b_0_ds;
g3_0 = x_b_0;

J0 = [g1_0, g2_0, g3_0];
g_inv_0 = inv(J0);

L3 = g3_0 / norm(g3_0);
L1_raw = cross(g2_0, L3);
if norm(L1_raw) < 1e-10
    L1_raw = cross([1; 0; 0], L3);
end
L1 = L1_raw / norm(L1_raw);
L2 = cross(L3, L1);

dx_m_t_dr = Xt' * dh_dr;
dx_m_t_ds = Xt' * dh_ds;

dx_b_t_dr = 0.5 * (Vnt .* thick)' * dh_dr;
dx_b_t_ds = 0.5 * (Vnt .* thick)' * dh_ds;

g1_t = dx_m_t_dr + zeta * dx_b_t_dr;
g2_t = dx_m_t_ds + zeta * dx_b_t_ds;
g3_t = 0.5 * (Vnt .* thick)' * h;

g_cov_0 = [dot(g1_0, g1_0); dot(g2_0, g2_0); dot(g1_0, g2_0); dot(g1_0, g3_0); dot(g2_0, g3_0)];
g_cov_t = [dot(g1_t, g1_t); dot(g2_t, g2_t); dot(g1_t, g2_t); dot(g1_t, g3_t); dot(g2_t, g3_t)];

B = zeros(5, 20);

for i = 1:4
    dof_idx = (i-1)*5 + (1:5);

    hi = h(i);
    dhi_dr = dh_dr(i);
    dhi_ds = dh_ds(i);
    a_i = thick(i);

    B(1, dof_idx(1:3)) = dhi_dr * g1_t';
    B(2, dof_idx(1:3)) = dhi_ds * g2_t';
    B(3, dof_idx(1:3)) = dhi_ds * g1_t' + dhi_dr * g2_t';
    B(4, dof_idx(1:3)) = dhi_dr * g3_t';
    B(5, dof_idx(1:3)) = dhi_ds * g3_t';

    V2_i = V2t(i, :)';
    V1_i = V1t(i, :)';

    vec_alpha_dr = 0.5 * a_i * dhi_dr * (-V2_i);
    vec_beta_dr  = 0.5 * a_i * dhi_dr * (V1_i);

    vec_alpha_ds = 0.5 * a_i * dhi_ds * (-V2_i);
    vec_beta_ds  = 0.5 * a_i * dhi_ds * (V1_i);

    vec_alpha_zeta = 0.5 * a_i * hi * (-V2_i);
    vec_beta_zeta  = 0.5 * a_i * hi * (V1_i);

    B(1, dof_idx(4)) = zeta * dot(g1_t, vec_alpha_dr);
    B(1, dof_idx(5)) = zeta * dot(g1_t, vec_beta_dr);

    B(2, dof_idx(4)) = zeta * dot(g2_t, vec_alpha_ds);
    B(2, dof_idx(5)) = zeta * dot(g2_t, vec_beta_ds);

    B(3, dof_idx(4)) = zeta * (dot(g1_t, vec_alpha_ds) + dot(g2_t, vec_alpha_dr));
    B(3, dof_idx(5)) = zeta * (dot(g1_t, vec_beta_ds)  + dot(g2_t, vec_beta_dr));

    B(4, dof_idx(4)) = dot(g1_t, vec_alpha_zeta) + zeta * dot(g3_t, vec_alpha_dr);
    B(4, dof_idx(5)) = dot(g1_t, vec_beta_zeta)  + zeta * dot(g3_t, vec_beta_dr);

    B(5, dof_idx(4)) = dot(g2_t, vec_alpha_zeta) + zeta * dot(g3_t, vec_alpha_ds);
    B(5, dof_idx(5)) = dot(g2_t, vec_beta_zeta)  + zeta * dot(g3_t, vec_beta_ds);
end

end


function T = compute_covariant_to_local_transform(L1, L2, L3, g_inv_0)
g1_inv = g_inv_0(1, :)';
g2_inv = g_inv_0(2, :)';
g3_inv = g_inv_0(3, :)';

L = [L1'; L2'; L3'];
G_inv = [g1_inv, g2_inv, g3_inv];
l = L * G_inv;

T = zeros(5, 5);

T(1,1) = l(1,1)^2;
T(1,2) = l(1,2)^2;
T(1,3) = l(1,1)*l(1,2);

T(2,1) = l(2,1)^2;
T(2,2) = l(2,2)^2;
T(2,3) = l(2,1)*l(2,2);

T(3,1) = 2*l(1,1)*l(2,1);
T(3,2) = 2*l(1,2)*l(2,2);
T(3,3) = l(1,1)*l(2,2) + l(1,2)*l(2,1);

T(4,4) = l(1,1)*l(3,3);
T(4,5) = l(1,2)*l(3,3);

T(5,4) = l(2,1)*l(3,3);
T(5,5) = l(2,2)*l(3,3);

end


function E_cov = compute_covariant_GL_strain(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, ...
    is_mitc4_plus, aA_0, aB_0, aC_0, aD_0, aE_0, aA_t, aB_t, aC_t, aD_t, aE_t, ...
    g_rt_tied_0, g_st_tied_0, g_rt_tied_t, g_st_tied_t, ...
    g_rr_mb_0, g_ss_mb_0, g_rs_mb_0, g_rr_mb_t, g_ss_mb_t, g_rs_mb_t)

[~, ~, g_cov_0, g_cov_t] = eval_kinematics(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick);

E_cov = 0.5 * (g_cov_t - g_cov_0);

E_rt_A = 0.5 * (g_rt_tied_t(1) - g_rt_tied_0(1));
E_rt_B = 0.5 * (g_rt_tied_t(2) - g_rt_tied_0(2));
E_st_C = 0.5 * (g_st_tied_t(3) - g_st_tied_0(3));
E_st_D = 0.5 * (g_st_tied_t(4) - g_st_tied_0(4));

E_cov(4) = 0.5*(1+s)*E_rt_A + 0.5*(1-s)*E_rt_B;
E_cov(5) = 0.5*(1+r)*E_st_C + 0.5*(1-r)*E_st_D;

if is_mitc4_plus
    E_rr_A = 0.5 * (g_rr_mb_t(1) - g_rr_mb_0(1));
    E_rr_B = 0.5 * (g_rr_mb_t(2) - g_rr_mb_0(2));
    E_ss_C = 0.5 * (g_ss_mb_t(3) - g_ss_mb_0(3));
    E_ss_D = 0.5 * (g_ss_mb_t(4) - g_ss_mb_0(4));
    E_rs_E = 0.5 * (g_rs_mb_t(5) - g_rs_mb_0(5));

    E_rr_m = 0.5*(1 - 2*aA_t + s + 2*aA_t*s^2)*E_rr_A ...
           + 0.5*(1 - 2*aB_t - s + 2*aB_t*s^2)*E_rr_B ...
           + aC_t*(-1 + s^2)*E_ss_C ...
           + aD_t*(-1 + s^2)*E_ss_D ...
           + aE_t*(-1 + s^2)*E_rs_E;

    E_ss_m = aA_t*(-1 + r^2)*E_rr_A ...
           + aB_t*(-1 + r^2)*E_rr_B ...
           + 0.5*(1 - 2*aC_t + r + 2*aC_t*r^2)*E_ss_C ...
           + 0.5*(1 - 2*aD_t - r + 2*aD_t*r^2)*E_ss_D ...
           + aE_t*(-1 + r^2)*E_rs_E;

    E_rs_m = 0.25*(r + 4*aA_t*r*s)*E_rr_A ...
           + 0.25*(-r + 4*aB_t*r*s)*E_rr_B ...
           + 0.25*(s + 4*aC_t*r*s)*E_ss_C ...
           + 0.25*(-s + 4*aD_t*r*s)*E_ss_D ...
           + (1 + aE_t*r*s)*E_rs_E;

    [~, ~, g_cov_0_m, g_cov_t_m] = eval_kinematics(r, s, 0, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick);
    E_cov_m_direct = 0.5 * (g_cov_t_m(1:3) - g_cov_0_m(1:3));

    E_cov(1:3) = E_cov(1:3) - E_cov_m_direct + [E_rr_m; E_ss_m; E_rs_m];
end

E_cov(3:5) = 2 * E_cov(3:5);

end


function K_geo = compute_geometric_stiffness(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, S_local, T, g_inv_0)
K_geo = zeros(20, 20);

S_cov = T' * S_local;
S_rr = S_cov(1);
S_ss = S_cov(2);
S_rs = S_cov(3);
S_rt = S_cov(4);
S_st = S_cov(5);

h = 0.25 * [(1-r)*(1-s); (1+r)*(1-s); (1+r)*(1+s); (1-r)*(1+s)];
dh_dr = 0.25 * [-(1-s);  (1-s);  (1+s); -(1+s)];
dh_ds = 0.25 * [-(1-r); -(1+r);  (1+r);  (1-r)];

for i = 1:4
    idx_i_trans = (i-1)*5 + (1:3);
    dhi_dr = dh_dr(i);
    dhi_ds = dh_ds(i);
    a_i = thick(i);

    for j = 1:4
        idx_j_trans = (j-1)*5 + (1:3);
        dhj_dr = dh_dr(j);
        dhj_ds = dh_ds(j);

        N_trans = S_rr * (dhi_dr * dhj_dr) ...
                + S_ss * (dhi_ds * dhj_ds) ...
                + S_rs * (dhi_dr * dhj_ds + dhi_ds * dhj_dr);

        K_geo(idx_i_trans, idx_j_trans) = K_geo(idx_i_trans, idx_j_trans) + N_trans * eye(3);

        if i == j
            N_rot_quad = 0.5 * a_i * (S_rt * (dhi_dr * h(i)) + S_st * (dhi_ds * h(i)));
            K_geo((i-1)*5 + 4, (i-1)*5 + 4) = K_geo((i-1)*5 + 4, (i-1)*5 + 4) + N_rot_quad;
            K_geo((i-1)*5 + 5, (i-1)*5 + 5) = K_geo((i-1)*5 + 5, (i-1)*5 + 5) + N_rot_quad;
        end
    end
end

end
