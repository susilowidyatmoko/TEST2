import numpy as np
from mitc4_element import mitc4_plus_element, material_matrix, update_geometry

def nonlinear_solver(X0, elements, thick, E, nu, fixed_dofs, load_dofs, total_load, num_steps=10, is_mitc4_plus=True, Vn0=None, max_iter=25, tol=1e-4):
    Nnodes = len(X0)
    Nelems = len(elements)
    Ndofs = 5 * Nnodes

    if np.isscalar(thick):
        thick = thick * np.ones(Nnodes)

    if Vn0 is None:
        Vn0 = np.tile([0.0, 0.0, 1.0], (Nnodes, 1))

    V10 = np.zeros((Nnodes, 3))
    V20 = np.zeros((Nnodes, 3))

    for i in range(Nnodes):
        Vn0[i] /= np.linalg.norm(Vn0[i])
        ref_v = np.array([1.0, 0.0, 0.0]) if abs(Vn0[i, 0]) < 0.8 else np.array([0.0, 1.0, 0.0])
        v1_tmp = np.cross(ref_v, Vn0[i])
        V10[i] = v1_tmp / np.linalg.norm(v1_tmp)
        v2_tmp = np.cross(Vn0[i], V10[i])
        V20[i] = v2_tmp / np.linalg.norm(v2_tmp)

    C_mat = material_matrix(E, nu)

    Xt = X0.copy()
    Vnt = Vn0.copy()
    V1t = V10.copy()
    V2t = V20.copy()
    U_total = np.zeros(Ndofs)

    all_dofs = np.arange(Ndofs)
    free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    history_disp = []
    history_loads = []

    for step in range(1, num_steps + 1):
        step_factor = step / float(num_steps)
        F_ext = np.zeros(Ndofs)
        F_ext[load_dofs] = step_factor * total_load

        iter_count = 0
        converged = False
        rel_res = 1.0

        while not converged and iter_count < max_iter:
            iter_count += 1

            K_glob = np.zeros((Ndofs, Ndofs))
            F_int = np.zeros(Ndofs)

            for e in range(Nelems):
                nodes_e = elements[e]
                dof_e = []
                for n in range(4):
                    dof_e.extend(range(nodes_e[n]*5, (nodes_e[n]+1)*5))
                dof_e = np.array(dof_e)

                X0_e = X0[nodes_e]
                Xt_e = Xt[nodes_e]
                Vn0_e = Vn0[nodes_e]; V10_e = V10[nodes_e]; V20_e = V20[nodes_e]
                Vnt_e = Vnt[nodes_e]; V1t_e = V1t[nodes_e]; V2t_e = V2t[nodes_e]
                thick_e = thick[nodes_e]

                Ke, Fe = mitc4_plus_element(X0_e, Xt_e, Vn0_e, V10_e, V20_e, Vnt_e, V1t_e, V2t_e, thick_e, C_mat, is_mitc4_plus)

                K_glob[np.ix_(dof_e, dof_e)] += Ke
                F_int[dof_e] += Fe

            Residual = F_ext - F_int
            res_norm = np.linalg.norm(Residual[free_dofs])
            ext_norm = np.linalg.norm(F_ext[free_dofs])
            if ext_norm < 1e-10:
                ext_norm = 1.0
            rel_res = res_norm / ext_norm

            if rel_res < tol and iter_count > 1:
                converged = True
                break

            dU_free = np.linalg.solve(K_glob[np.ix_(free_dofs, free_dofs)], Residual[free_dofs])
            dU = np.zeros(Ndofs)
            dU[free_dofs] = dU_free

            Xt, Vnt, V1t, V2t = update_geometry(Xt, Vnt, V1t, V2t, dU)
            U_total += dU

        print(f"Step {step:2d}/{num_steps:2d}: Iterations = {iter_count:2d}, Rel Residual = {rel_res:.4e}")
        history_disp.append(U_total.copy())
        history_loads.append(F_ext[load_dofs].copy())

    return np.array(history_disp), np.array(history_loads)
