import numpy as np

def material_matrix(E, nu, k_shear=5.0/6.0):
    factor = E / (1.0 - nu**2)
    G = E / (2.0 * (1.0 + nu))
    C = np.zeros((5, 5))
    C[0, 0] = factor
    C[0, 1] = factor * nu
    C[1, 0] = factor * nu
    C[1, 1] = factor
    C[2, 2] = factor * (1.0 - nu) / 2.0
    C[3, 3] = k_shear * G
    C[4, 4] = k_shear * G
    return C

def update_geometry(X, V_n, V_1, V_2, U):
    Nnodes = len(X)
    X_new = X.copy()
    V_n_new = V_n.copy()
    V_1_new = V_1.copy()
    V_2_new = V_2.copy()

    for i in range(Nnodes):
        idx = slice(i * 5, (i + 1) * 5)
        u_i = U[idx][:3]
        alpha_i = U[idx][3]
        beta_i = U[idx][4]

        X_new[i] += u_i

        theta_vec = -alpha_i * V_2[i] + beta_i * V_1[i]
        theta_mag = np.linalg.norm(theta_vec)

        if theta_mag > 1e-12:
            axis_rot = theta_vec / theta_mag
            q0 = np.cos(theta_mag / 2.0)
            q1 = axis_rot[0] * np.sin(theta_mag / 2.0)
            q2 = axis_rot[1] * np.sin(theta_mag / 2.0)
            q3 = axis_rot[2] * np.sin(theta_mag / 2.0)

            Q = 2.0 * np.array([
                [q0**2 + q1**2 - 0.5, q1*q2 - q0*q3, q3*q1 + q0*q2],
                [q1*q2 + q0*q3, q0**2 + q2**2 - 0.5, q2*q3 - q0*q1],
                [q3*q1 - q0*q2, q2*q3 + q0*q1, q0**2 + q3**2 - 0.5]
            ])

            V_n_new[i] = Q @ V_n[i]
            V_1_new[i] = Q @ V_1[i]
            V_2_new[i] = Q @ V_2[i]

        V_n_new[i] /= np.linalg.norm(V_n_new[i])

        V_1_proj = V_1_new[i] - np.dot(V_1_new[i], V_n_new[i]) * V_n_new[i]
        if np.linalg.norm(V_1_proj) > 1e-8:
            V_1_new[i] = V_1_proj / np.linalg.norm(V_1_proj)
        else:
            ref_v = np.array([1.0, 0.0, 0.0]) if abs(V_n_new[i, 0]) < 0.8 else np.array([0.0, 1.0, 0.0])
            v1_tmp = np.cross(ref_v, V_n_new[i])
            V_1_new[i] = v1_tmp / np.linalg.norm(v1_tmp)

        v2_tmp = np.cross(V_n_new[i], V_1_new[i])
        V_2_new[i] = v2_tmp / np.linalg.norm(v2_tmp)

    return X_new, V_n_new, V_1_new, V_2_new

def eval_kinematics(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick):
    h = 0.25 * np.array([(1-r)*(1-s), (1+r)*(1-s), (1+r)*(1+s), (1-r)*(1+s)])
    dh_dr = 0.25 * np.array([-(1-s), (1-s), (1+s), -(1+s)])
    dh_ds = 0.25 * np.array([-(1-r), -(1+r), (1+r), (1-r)])

    dx_m_0_dr = X0.T @ dh_dr
    dx_m_0_ds = X0.T @ dh_ds

    x_b_0 = 0.5 * (Vn0 * thick[:, None]).T @ h
    dx_b_0_dr = 0.5 * (Vn0 * thick[:, None]).T @ dh_dr
    dx_b_0_ds = 0.5 * (Vn0 * thick[:, None]).T @ dh_ds

    g1_0 = dx_m_0_dr + zeta * dx_b_0_dr
    g2_0 = dx_m_0_ds + zeta * dx_b_0_ds
    g3_0 = x_b_0

    J0 = np.column_stack([g1_0, g2_0, g3_0])
    g_inv_0 = np.linalg.inv(J0)

    L3 = g3_0 / np.linalg.norm(g3_0)
    L1_raw = np.cross(g2_0, L3)
    if np.linalg.norm(L1_raw) < 1e-10:
        L1_raw = np.cross(np.array([1.0, 0.0, 0.0]), L3)
    L1 = L1_raw / np.linalg.norm(L1_raw)
    L2 = np.cross(L3, L1)

    dx_m_t_dr = Xt.T @ dh_dr
    dx_m_t_ds = Xt.T @ dh_ds
    dx_b_t_dr = 0.5 * (Vnt * thick[:, None]).T @ dh_dr
    dx_b_t_ds = 0.5 * (Vnt * thick[:, None]).T @ dh_ds

    g1_t = dx_m_t_dr + zeta * dx_b_t_dr
    g2_t = dx_m_t_ds + zeta * dx_b_t_ds
    g3_t = 0.5 * (Vnt * thick[:, None]).T @ h

    g_cov_0 = np.array([np.dot(g1_0, g1_0), np.dot(g2_0, g2_0), np.dot(g1_0, g2_0), np.dot(g1_0, g3_0), np.dot(g2_0, g3_0)])
    g_cov_t = np.array([np.dot(g1_t, g1_t), np.dot(g2_t, g2_t), np.dot(g1_t, g2_t), np.dot(g1_t, g3_t), np.dot(g2_t, g3_t)])

    B = np.zeros((5, 20))
    for i in range(4):
        dof_idx = slice(i * 5, (i + 1) * 5)
        hi = h[i]
        dhi_dr = dh_dr[i]
        dhi_ds = dh_ds[i]
        a_i = thick[i]

        B[0, dof_idx][:3] = dhi_dr * g1_t
        B[1, dof_idx][:3] = dhi_ds * g2_t
        B[2, dof_idx][:3] = dhi_ds * g1_t + dhi_dr * g2_t
        B[3, dof_idx][:3] = dhi_dr * g3_t
        B[4, dof_idx][:3] = dhi_ds * g3_t

        V2_i = V2t[i]
        V1_i = V1t[i]

        vec_alpha_dr = 0.5 * a_i * dhi_dr * (-V2_i)
        vec_beta_dr  = 0.5 * a_i * dhi_dr * (V1_i)
        vec_alpha_ds = 0.5 * a_i * dhi_ds * (-V2_i)
        vec_beta_ds  = 0.5 * a_i * dhi_ds * (V1_i)
        vec_alpha_zeta = 0.5 * a_i * hi * (-V2_i)
        vec_beta_zeta  = 0.5 * a_i * hi * (V1_i)

        B[0, dof_idx][3] = zeta * np.dot(g1_t, vec_alpha_dr)
        B[0, dof_idx][4] = zeta * np.dot(g1_t, vec_beta_dr)

        B[1, dof_idx][3] = zeta * np.dot(g2_t, vec_alpha_ds)
        B[1, dof_idx][4] = zeta * np.dot(g2_t, vec_beta_ds)

        B[2, dof_idx][3] = zeta * (np.dot(g1_t, vec_alpha_ds) + np.dot(g2_t, vec_alpha_dr))
        B[2, dof_idx][4] = zeta * (np.dot(g1_t, vec_beta_ds) + np.dot(g2_t, vec_beta_dr))

        B[3, dof_idx][3] = np.dot(g1_t, vec_alpha_zeta) + zeta * np.dot(g3_t, vec_alpha_dr)
        B[3, dof_idx][4] = np.dot(g1_t, vec_beta_zeta) + zeta * np.dot(g3_t, vec_beta_dr)

        B[4, dof_idx][3] = np.dot(g2_t, vec_alpha_zeta) + zeta * np.dot(g3_t, vec_alpha_ds)
        B[4, dof_idx][4] = np.dot(g2_t, vec_beta_zeta) + zeta * np.dot(g3_t, vec_beta_ds)

    return J0, B, g_cov_0, g_cov_t, g_inv_0, L1, L2, L3

def compute_covariant_to_local_transform(L1, L2, L3, g_inv_0):
    l = np.array([L1, L2, L3]) @ g_inv_0.T
    T = np.zeros((5, 5))
    T[0, 0] = l[0, 0]**2
    T[0, 1] = l[0, 1]**2
    T[0, 2] = l[0, 0] * l[0, 1]
    T[1, 0] = l[1, 0]**2
    T[1, 1] = l[1, 1]**2
    T[1, 2] = l[1, 0] * l[1, 1]
    T[2, 0] = 2.0 * l[0, 0] * l[1, 0]
    T[2, 1] = 2.0 * l[0, 1] * l[1, 1]
    T[2, 2] = l[0, 0] * l[1, 1] + l[0, 1] * l[1, 0]
    T[3, 3] = l[0, 0] * l[2, 2]
    T[3, 4] = l[0, 1] * l[2, 2]
    T[4, 3] = l[1, 0] * l[2, 2]
    T[4, 4] = l[1, 1] * l[2, 2]
    return T

def compute_geometric_stiffness(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, S_local, T, g_inv_0):
    K_geo = np.zeros((20, 20))
    S_cov = T.T @ S_local
    S_rr, S_ss, S_rs, S_rt, S_st = S_cov[0], S_cov[1], S_cov[2], S_cov[3], S_cov[4]

    h = 0.25 * np.array([(1-r)*(1-s), (1+r)*(1-s), (1+r)*(1+s), (1-r)*(1+s)])
    dh_dr = 0.25 * np.array([-(1-s), (1-s), (1+s), -(1+s)])
    dh_ds = 0.25 * np.array([-(1-r), -(1+r), (1+r), (1-r)])

    for i in range(4):
        trans_i = slice(i * 5, i * 5 + 3)
        dhi_dr, dhi_ds = dh_dr[i], dh_ds[i]
        a_i = thick[i]

        for j in range(4):
            trans_j = slice(j * 5, j * 5 + 3)
            dhj_dr, dhj_ds = dh_dr[j], dh_ds[j]

            N_trans = S_rr * (dhi_dr * dhj_dr) + S_ss * (dhi_ds * dhj_ds) + S_rs * (dhi_dr * dhj_ds + dhi_ds * dhj_dr)
            K_geo[trans_i, trans_j] += N_trans * np.eye(3)

            if i == j:
                N_rot = 0.5 * a_i * (S_rt * (dhi_dr * h[i]) + S_st * (dhi_ds * h[i]))
                K_geo[i * 5 + 3, i * 5 + 3] += N_rot
                K_geo[i * 5 + 4, i * 5 + 4] += N_rot

    return K_geo

def mitc4_plus_element(X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, C_mat, is_mitc4_plus=True):
    Ke = np.zeros((20, 20))
    Fe = np.zeros(20)

    if np.isscalar(thick):
        thick = thick * np.ones(4)

    gp = np.array([-1.0/np.sqrt(3), 1.0/np.sqrt(3)])
    gw = np.array([1.0, 1.0])

    r_tying_shear = [0, 0, -1, 1]
    s_tying_shear = [-1, 1, 0, 0]

    B_shear_tied = []

    for i in range(4):
        J0, B_sh, g_cov_0, g_cov_t, g_inv_0, L1, L2, L3 = eval_kinematics(
            r_tying_shear[i], s_tying_shear[i], 0.0, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick
        )
        B_shear_tied.append(B_sh)

    for ig in range(2):
        r = gp[ig]
        wr = gw[ig]
        for jg in range(2):
            s = gp[jg]
            ws = gw[jg]
            for kg in range(2):
                zeta = gp[kg]
                wz = gw[kg]

                J0, B_direct, g_cov_0, g_cov_t, g_inv_0, L1, L2, L3 = eval_kinematics(
                    r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick
                )

                detJ0 = np.linalg.det(J0)
                dV = detJ0 * wr * ws * wz

                B_cov = B_direct.copy()
                B_cov[3, :] = 0.5 * (1.0 + s) * B_shear_tied[0][3, :] + 0.5 * (1.0 - s) * B_shear_tied[1][3, :]
                B_cov[4, :] = 0.5 * (1.0 + r) * B_shear_tied[2][4, :] + 0.5 * (1.0 - r) * B_shear_tied[3][4, :]

                T = compute_covariant_to_local_transform(L1, L2, L3, g_inv_0)
                B_local = T @ B_cov

                E_cov = 0.5 * (g_cov_t - g_cov_0)
                E_cov[2:] *= 2.0
                E_local = T @ E_cov

                S_vec = C_mat @ E_local

                K_geo = compute_geometric_stiffness(r, s, zeta, X0, Xt, Vn0, V10, V20, Vnt, V1t, V2t, thick, S_vec, T, g_inv_0)

                Ke += (B_local.T @ C_mat @ B_local + K_geo) * dV
                Fe += (B_local.T @ S_vec) * dV

    return Ke, Fe
