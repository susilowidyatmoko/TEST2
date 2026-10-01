import numpy as np

def mesh_cantilever(L, b, nx, ny):
    x = np.linspace(0, L, nx + 1)
    y = np.linspace(-b/2, b/2, ny + 1)
    YY, XX = np.meshgrid(y, x)
    X = np.column_stack([XX.ravel(), YY.ravel(), np.zeros(XX.size)])

    elements = []
    for i in range(nx):
        for j in range(ny):
            n1 = i * (ny + 1) + j
            n2 = (i + 1) * (ny + 1) + j
            n3 = (i + 1) * (ny + 1) + j + 1
            n4 = i * (ny + 1) + j + 1
            elements.append([n1, n2, n3, n4])

    Vn0 = np.tile([0.0, 0.0, 1.0], (len(X), 1))
    return X, np.array(elements, dtype=int), Vn0

def mesh_slit_annular(Ri, Ro, nr, ntheta):
    r = np.linspace(Ri, Ro, nr + 1)
    theta = np.linspace(0, 2 * np.pi, ntheta + 1)
    TT, RR = np.meshgrid(theta, r)
    XX = RR * np.cos(TT)
    YY = RR * np.sin(TT)
    ZZ = np.zeros_like(XX)
    X = np.column_stack([XX.ravel(), YY.ravel(), ZZ.ravel()])

    elements = []
    for j in range(ntheta):
        for i in range(nr):
            n1 = j * (nr + 1) + i
            n2 = j * (nr + 1) + i + 1
            n3 = (j + 1) * (nr + 1) + i + 1
            n4 = (j + 1) * (nr + 1) + i
            elements.append([n1, n2, n3, n4])

    Vn0 = np.tile([0.0, 0.0, 1.0], (len(X), 1))
    return X, np.array(elements, dtype=int), Vn0

def mesh_hemisphere(R, phi_hole_deg, n_phi, n_theta):
    phi_min = np.radians(phi_hole_deg)
    phi_max = np.pi / 2
    theta_min = 0
    theta_max = np.pi / 2

    phi = np.linspace(phi_min, phi_max, n_phi + 1)
    theta = np.linspace(theta_min, theta_max, n_theta + 1)
    TT, PP = np.meshgrid(theta, phi)

    XX = R * np.sin(PP) * np.cos(TT)
    YY = R * np.sin(PP) * np.sin(TT)
    ZZ = R * np.cos(PP)
    X = np.column_stack([XX.ravel(), YY.ravel(), ZZ.ravel()])

    elements = []
    for j in range(n_theta):
        for i in range(n_phi):
            n1 = j * (n_phi + 1) + i
            n2 = j * (n_phi + 1) + i + 1
            n3 = (j + 1) * (n_phi + 1) + i + 1
            n4 = (j + 1) * (n_phi + 1) + i
            elements.append([n1, n2, n3, n4])

    norms = np.linalg.norm(X, axis=1, keepdims=True)
    Vn0 = np.where(norms > 1e-12, X / norms, np.array([[0.0, 0.0, 1.0]]))
    return X, np.array(elements, dtype=int), Vn0

def mesh_scordelis(R, L, theta_deg, nx, ny):
    theta_rad = np.radians(theta_deg)
    x = np.linspace(0, L / 2, nx + 1)
    theta = np.linspace(0, theta_rad, ny + 1)
    TT, XX = np.meshgrid(theta, x)
    YY = R * np.sin(TT)
    ZZ = R * np.cos(TT)
    X = np.column_stack([XX.ravel(), YY.ravel(), ZZ.ravel()])

    elements = []
    for j in range(ny):
        for i in range(nx):
            n1 = j * (nx + 1) + i
            n2 = j * (nx + 1) + i + 1
            n3 = (j + 1) * (nx + 1) + i + 1
            n4 = (j + 1) * (nx + 1) + i
            elements.append([n1, n2, n3, n4])

    TT_vec = TT.ravel()
    Vn0 = np.column_stack([np.zeros_like(TT_vec), np.sin(TT_vec), np.cos(TT_vec)])
    return X, np.array(elements, dtype=int), Vn0

def mesh_cylinder(R, L, theta_deg, nx, ny, is_distorted=False):
    theta_rad = np.radians(theta_deg)
    if not is_distorted:
        x = np.linspace(0, L, nx + 1)
        s = np.linspace(0, theta_rad, ny + 1)
    else:
        rx = np.arange(1, nx + 1)
        x = L * np.concatenate([[0], np.cumsum(rx) / np.sum(rx)])
        rs = np.arange(1, ny + 1)
        s = theta_rad * np.concatenate([[0], np.cumsum(rs) / np.sum(rs)])

    SS, XX = np.meshgrid(s, x)
    YY = R * np.sin(SS)
    ZZ = R * np.cos(SS)
    X = np.column_stack([XX.ravel(), YY.ravel(), ZZ.ravel()])

    elements = []
    for j in range(ny):
        for i in range(nx):
            n1 = j * (nx + 1) + i
            n2 = j * (nx + 1) + i + 1
            n3 = (j + 1) * (nx + 1) + i + 1
            n4 = (j + 1) * (nx + 1) + i
            elements.append([n1, n2, n3, n4])

    SS_vec = SS.ravel()
    Vn0 = np.column_stack([np.zeros_like(SS_vec), np.sin(SS_vec), np.cos(SS_vec)])
    return X, np.array(elements, dtype=int), Vn0

def mesh_twisted_cantilever(L, w, total_twist_deg, nx, ny):
    x = np.linspace(0, L, nx + 1)
    y = np.linspace(-w / 2, w / 2, ny + 1)
    YY, XX = np.meshgrid(y, x)
    theta_x = np.radians(total_twist_deg) * (XX / L)

    XX_t = XX
    YY_t = YY * np.cos(theta_x)
    ZZ_t = YY * np.sin(theta_x)
    X = np.column_stack([XX_t.ravel(), YY_t.ravel(), ZZ_t.ravel()])

    elements = []
    for j in range(ny):
        for i in range(nx):
            n1 = j * (nx + 1) + i
            n2 = j * (nx + 1) + i + 1
            n3 = (j + 1) * (nx + 1) + i + 1
            n4 = (j + 1) * (nx + 1) + i
            elements.append([n1, n2, n3, n4])

    ang_vec = theta_x.ravel()
    Vn0 = np.column_stack([np.zeros_like(ang_vec), -np.sin(ang_vec), np.cos(ang_vec)])
    return X, np.array(elements, dtype=int), Vn0

def mesh_hook(width, s_total, n_s, n_w):
    R1 = 14.12
    phi1 = np.radians(150)
    R2 = 4.6
    phi2 = np.radians(180)

    L1 = R1 * phi1
    L2 = R2 * phi2
    L_tot = L1 + L2

    s = np.linspace(0, L_tot, n_s + 1)
    z = np.linspace(-width / 2, width / 2, n_w + 1)

    X_pts = np.zeros((n_s + 1, 3))
    Vn_pts = np.zeros((n_s + 1, 3))
    for i, si in enumerate(s):
        if si <= L1:
            ang = si / R1
            X_pts[i, 0] = R1 * (1 - np.cos(ang))
            X_pts[i, 1] = R1 * np.sin(ang)
            Vn_pts[i] = [-np.cos(ang), np.sin(ang), 0.0]
        else:
            s2 = si - L1
            ang2 = s2 / R2
            C1 = [R1 * (1 - np.cos(phi1)), R1 * np.sin(phi1)]
            X_pts[i, 0] = C1[0] + R2 * (np.sin(phi1 - ang2) - np.sin(phi1))
            X_pts[i, 1] = C1[1] - R2 * (np.cos(phi1 - ang2) - np.cos(phi1))
            Vn_pts[i] = [np.cos(phi1 - ang2), np.sin(phi1 - ang2), 0.0]

    ZZ, XX_line = np.meshgrid(z, X_pts[:, 0])
    _, YY_line = np.meshgrid(z, X_pts[:, 1])
    _, VnX_line = np.meshgrid(z, Vn_pts[:, 0])
    _, VnY_line = np.meshgrid(z, Vn_pts[:, 1])
    _, VnZ_line = np.meshgrid(z, Vn_pts[:, 2])

    X = np.column_stack([XX_line.ravel(), YY_line.ravel(), ZZ.ravel()])
    Vn0 = np.column_stack([VnX_line.ravel(), VnY_line.ravel(), VnZ_line.ravel()])

    elements = []
    for j in range(n_w):
        for i in range(n_s):
            n1 = j * (n_s + 1) + i
            n2 = j * (n_s + 1) + i + 1
            n3 = (j + 1) * (n_s + 1) + i + 1
            n4 = (j + 1) * (n_s + 1) + i
            elements.append([n1, n2, n3, n4])
    return X, np.array(elements, dtype=int), Vn0
