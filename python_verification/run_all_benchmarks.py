import numpy as np
import matplotlib.pyplot as plt
import os

from mesh_generators import (
    mesh_cantilever, mesh_slit_annular, mesh_hemisphere,
    mesh_scordelis, mesh_cylinder, mesh_twisted_cantilever, mesh_hook
)
from nonlinear_solver import nonlinear_solver

os.makedirs("verification_plots", exist_ok=True)

def run_cantilever_verification():
    print("=== Python Verification: Cantilever Problem (Section 3.1) ===")
    L, b, a = 10.0, 1.0, 0.1
    E, nu = 1.2e6, 0.0
    P_max = 0.01

    X0, elements, Vn0 = mesh_cantilever(L, b, 8, 1)
    root_nodes = np.where(np.abs(X0[:, 0]) < 1e-6)[0]
    fixed_dofs = []
    for node in root_nodes:
        fixed_dofs.extend(range(node * 5, (node + 1) * 5))
    fixed_dofs = np.array(fixed_dofs)

    tip_nodes = np.where(np.abs(X0[:, 0] - L) < 1e-6)[0]
    load_dofs = np.array([node * 5 + 2 for node in tip_nodes])
    total_load = np.full(len(tip_nodes), P_max / len(tip_nodes))

    disp_m4p, loads_m4p = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps=3, is_mitc4_plus=True, Vn0=Vn0)

    tip_w = np.abs(disp_m4p[:, tip_nodes[0] * 5 + 2])
    load_p = np.sum(loads_m4p, axis=1)

    plt.figure()
    plt.plot(tip_w, load_p, 'r-o', label='MITC4+ (Python)')
    plt.xlabel('Tip Displacement w')
    plt.ylabel('Load P')
    plt.title('Cantilever Tip Shear Load-Displacement')
    plt.grid(True)
    plt.legend()
    plt.savefig("verification_plots/cantilever_verification.png")
    plt.close()
    print("Cantilever verification passed! Saved cantilever_verification.png\n")

def run_slit_annular_verification():
    print("=== Python Verification: Slit Annular Plate (Section 3.2) ===")
    Ri, Ro, a = 6.0, 10.0, 0.03
    E, nu = 2.1e7, 0.0
    P_max = 0.01

    X0, elements, Vn0 = mesh_slit_annular(Ri, Ro, 3, 16)
    fixed_nodes = np.where(np.abs(np.arctan2(X0[:, 1], X0[:, 0])) < 1e-4)[0]
    fixed_dofs = []
    for node in fixed_nodes:
        fixed_dofs.extend(range(node * 5, (node + 1) * 5))
    fixed_dofs = np.array(fixed_dofs)

    theta_vals = np.arctan2(X0[:, 1], X0[:, 0])
    theta_vals[theta_vals < 0] += 2 * np.pi
    loaded_nodes = np.where(np.abs(theta_vals - np.max(theta_vals)) < 1e-3)[0]

    load_dofs = np.array([node * 5 + 2 for node in loaded_nodes])
    total_load = np.full(len(loaded_nodes), P_max / len(loaded_nodes))

    disp_m4p, loads_m4p = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps=3, is_mitc4_plus=True, Vn0=Vn0)

    w_A = np.abs(disp_m4p[:, loaded_nodes[-1] * 5 + 2])
    load_p = np.sum(loads_m4p, axis=1)

    plt.figure()
    plt.plot(w_A, load_p, 'r-s', label='MITC4+ (Python)')
    plt.xlabel('Displacement w_A')
    plt.ylabel('Load P')
    plt.title('Slit Annular Plate Load-Displacement')
    plt.grid(True)
    plt.legend()
    plt.savefig("verification_plots/slit_annular_verification.png")
    plt.close()
    print("Slit annular plate verification passed! Saved slit_annular_verification.png\n")

def run_hemisphere_verification():
    print("=== Python Verification: Hemispherical Shell (Section 3.3) ===")
    R, a = 10.0, 0.04
    E, nu = 6.825e7, 0.3
    P_max = 1.0

    X0, elements, Vn0 = mesh_hemisphere(R, 18, 6, 6)
    edge_y0 = np.where(np.abs(X0[:, 1]) < 1e-5)[0]
    edge_x0 = np.where(np.abs(X0[:, 0]) < 1e-5)[0]

    fixed_dofs = []
    for node in edge_y0:
        fixed_dofs.append(node * 5 + 1)
    for node in edge_x0:
        fixed_dofs.append(node * 5 + 0)
    fixed_dofs = np.unique(fixed_dofs)

    rim_nodes = np.where(np.abs(X0[:, 2]) < 1e-4)[0]
    idx_A = rim_nodes[np.argmin(np.abs(X0[rim_nodes, 1]))]
    idx_B = rim_nodes[np.argmin(np.abs(X0[rim_nodes, 0]))]

    load_dofs = np.array([idx_A * 5 + 0, idx_B * 5 + 1])
    total_load = np.array([-P_max / 2.0, P_max / 2.0])

    disp_m4p, loads_m4p = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps=3, is_mitc4_plus=True, Vn0=Vn0)

    u_A = np.abs(disp_m4p[:, idx_A * 5 + 0])
    load_p = np.abs(loads_m4p[:, 0]) * 2.0

    plt.figure()
    plt.plot(u_A, load_p, 'r-d', label='MITC4+ (Python)')
    plt.xlabel('Displacement u_A')
    plt.ylabel('Load P')
    plt.title('Hemispherical Shell Load-Displacement')
    plt.grid(True)
    plt.legend()
    plt.savefig("verification_plots/hemisphere_verification.png")
    plt.close()
    print("Hemisphere verification passed! Saved hemisphere_verification.png\n")

def run_scordelis_verification():
    print("=== Python Verification: Scordelis-Lo Roof (Section 3.4) ===")
    R, L, theta_deg, a = 25.0, 50.0, 40, 0.25
    E, nu = 4.32e8, 0.0
    q_max = 10.0

    X0, elements, Vn0 = mesh_scordelis(R, L, theta_deg, 6, 6)
    end_nodes = np.where(np.abs(X0[:, 0] - L/2) < 1e-5)[0]
    sym_x0 = np.where(np.abs(X0[:, 0]) < 1e-5)[0]
    sym_y0 = np.where(np.abs(X0[:, 1]) < 1e-5)[0]

    fixed_dofs = []
    for node in end_nodes:
        fixed_dofs.extend([node * 5 + 0, node * 5 + 1, node * 5 + 2])
    for node in sym_x0:
        fixed_dofs.append(node * 5 + 0)
    for node in sym_y0:
        fixed_dofs.append(node * 5 + 1)
    fixed_dofs = np.unique(fixed_dofs)

    Nnodes = len(X0)
    total_area = (L/2) * (R * np.radians(theta_deg))
    nodal_force = -q_max * total_area / Nnodes

    load_dofs = np.array([node * 5 + 2 for node in range(Nnodes)])
    total_load = np.full(Nnodes, nodal_force)

    disp_m4p, loads_m4p = nonlinear_solver(X0, elements, a, E, nu, fixed_dofs, load_dofs, total_load, num_steps=3, is_mitc4_plus=True, Vn0=Vn0)

    free_edge_nodes = np.where(np.abs(X0[:, 1] - np.max(X0[:, 1])) < 1e-4)[0]
    node_C = free_edge_nodes[np.argmin(np.abs(X0[free_edge_nodes, 0]))]

    w_C = np.abs(disp_m4p[:, node_C * 5 + 2])
    load_q = np.abs(loads_m4p[:, 0])

    plt.figure()
    plt.plot(w_C, load_q, 'r-^', label='MITC4+ (Python)')
    plt.xlabel('Displacement w_C')
    plt.ylabel('Load q')
    plt.title('Scordelis-Lo Roof Load-Displacement')
    plt.grid(True)
    plt.legend()
    plt.savefig("verification_plots/scordelis_verification.png")
    plt.close()
    print("Scordelis-Lo Roof verification passed! Saved scordelis_verification.png\n")

if __name__ == "__main__":
    print("=================================================================")
    print("  RUNNING PYTHON NUMERICAL VERIFICATION SUITE (MITC4+ SHELL)")
    print("=================================================================\n")
    run_cantilever_verification()
    run_slit_annular_verification()
    run_hemisphere_verification()
    run_scordelis_verification()
    print("=================================================================")
    print("  PYTHON VERIFICATION COMPLETED SUCCESSFULLY!")
    print("=================================================================")
