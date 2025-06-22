import os
import numpy as np
from meshpy.triangle import MeshInfo, build


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import inspect, re
import sympy as sp
from sympy.printing.octave import octave_code


u, v = sp.symbols('u v')
_sympy_namespace = {name: obj for name, obj in sp.__dict__.items()}

def make_flag_func(umin, umax, vmin, vmax, tolu=0.05, tolv=0):
    def flag(u,v):
        if u < umin + tolu or u > umax - tolu or \
           v < vmin + tolv or v > vmax - tolv:
            return -1
        return -2
    return flag


def create_input_files(sim_name,
                       output_dir,
                       domain_type,
                       domain_params,
                       max_area=0.001,
                       tols=(0,0)):
    """
    Creates:
      {sim_name}_p-distmesh.dat   raw distmesh points
      {sim_name}_t-distmesh.dat   raw distmesh triangles
      {sim_name}_Vertices         MATLAB-style 4×N array
      {sim_name}_Faces            MATLAB-style 10×M array
    """
    os.makedirs(output_dir, exist_ok=True)
    # 1) build boundary
    if domain_type=='rectangle':
        umin,umax,vmin,vmax = domain_params
        pts    = [(umin,vmin),(umax,vmin),(umax,vmax),(umin,vmax)]
        facets = [(i,(i+1)%4) for i in range(4)]
    elif domain_type=='circle':
        cx,cy,r = domain_params
        angles  = np.linspace(0,2*np.pi,32,endpoint=False)
        pts     = [(cx + r*np.cos(a), cy + r*np.sin(a)) for a in angles]
        facets  = [(i,(i+1)%len(pts)) for i in range(len(pts))]
    elif domain_type=='polygon':
        pts     = list(domain_params)
        facets  = [(i,(i+1)%len(pts)) for i in range(len(pts))]
    else:
        raise ValueError(domain_type)

    # 2) refine
    info = MeshInfo()
    info.set_points(pts)
    info.set_facets(facets)
    mesh = build(info, refinement_func=lambda verts, area: area>max_area)
    P_raw = np.array(mesh.points)      # (N,2)
    T_raw = np.array(mesh.elements)    # (M,3)

    # 3) save distmesh files
    pfile = os.path.join(output_dir, f"{sim_name}_p-distmesh.dat")
    tfile = os.path.join(output_dir, f"{sim_name}_t-distmesh.dat")
    with open(pfile,'w') as f:
        # f.write(f"{len(P_raw)}\n")
        for (u,v) in P_raw:
            f.write(f"{u:.7f}\t{v:.7f}\n")
    with open(tfile,'w') as f:
        for (v0,v1,v2) in T_raw:
            f.write(f"{v0+1}\t{v1+1}\t{v2+1}\n")

    N = len(P_raw)
    V = np.zeros((4, N), dtype=float)
    V[0, :] = np.arange(N)
    V[1:3, :] = P_raw.T
    # apply your flag_func to each (u,v)
    tolu, tolv = tols
    flag_func    = make_flag_func(umin, umax, vmin, vmax, tolu=tolu, tolv=tolv)

    V[3, :] = [flag_func(u, v) for u, v in P_raw]

    # write out the Vertices file
    vfile = os.path.join(output_dir, f"{sim_name}_Vertices")
    with open(vfile, 'w') as f:
        f.write(f"{N}\n")
        for j in range(N):
            # MATLAB fprintf('%5.0f\t%10.7f\t%10.7f\t%5.0f\n', P(:,j))
            f.write(f"{int(V[0,j])}\t{V[1,j]:.7f}\t{V[2,j]:.7f}\t{int(V[3,j])}\n")
    # 4) build MATLAB-style P,T
    M = len(T_raw)

    # 4a) build an edge→triangles map
    edge_map = {}
    for tri_idx, tri in enumerate(T_raw):
        for k in range(3):
            v1 = tri[(k+1)%3]
            v2 = tri[(k+2)%3]
            key = tuple(sorted((v1,v2)))
            edge_map.setdefault(key, []).append(tri_idx)

    # 4b) initialize neighbor arrays
    neigh_tri  = -np.ones((M,3), dtype=int)
    neigh_vert = -np.ones((M,3), dtype=int)

    # 4c) fill them
    for i, tri in enumerate(T_raw):
        for k in range(3):
            # edge opposite vertex k in triangle i
            a, b = tri[(k+1)%3], tri[(k+2)%3]
            adj = edge_map.get(tuple(sorted((a,b))), [])
            # find the “other” triangle in adj
            other = next((j for j in adj if j != i), None)
            if other is not None:
                neigh_tri[i, k] = other
                # the vertex in 'other' that is *not* a or b
                opp = [v for v in T_raw[other] if v not in (a,b)][0]
                neigh_vert[i, k] = opp

    # 4d) stack into 10×M array
    #   row 0: triangle index
    #   rows 1–3: vertices
    #   rows 4–6: opposite-vertex indices from each neighbor
    #   rows 7–9: neighbor-triangle indices
    T = np.vstack([
        np.arange(M),           # row 0
        T_raw.T,                # rows 1–3
        neigh_vert.T,           # rows 4–6
        neigh_tri.T             # rows 7–9
    ]).astype(int)

    # 5) write out the Faces file
    ffile = os.path.join(output_dir, f"{sim_name}_Faces")
    with open(ffile, 'w') as f:
        f.write(f"{M}\n")
        for j in range(M):
            f.write("\t".join(str(x) for x in T[:,j]) + "\n")
    
    # 6) **Plot the triangle mesh** in the UV‐plane (inline, no saving)

    fig, ax = plt.subplots()
    ax.triplot(P_raw[:,0], P_raw[:,1], T_raw, color='black', lw=0.5)
    # ax.set_aspect('equal')
    ax.set_xlabel('u')
    ax.set_ylabel('v')
    ax.set_title(f'{sim_name} UV-mesh')
    plt.tight_layout()
    plt.show()
    


def make_in_file(params, output_dir, sim_name):
    lines = []
    lines.append(params['vertices_file'])
    lines.append(params['faces_file'])
    lines.append(str(params['loops']))
    lines.append(str(params['save_every']))

    # flatten abar, bbar and pos0, casting to strings
    lines.extend(str(x) for x in params['abar'])
    lines.extend(str(x) for x in params['bbar'])
    lines.extend([ str(params['thickness']),
                   str(params['E']),
                   str(params['nu']) ])
    lines.extend(str(x) for x in params['pos0'])
    lines.extend([ str(params['lambdaG']),
                   str(params['muG']) ])

    # γ‐table
    for k in (0,1):
        for i in (0,1):
            for j in (0,1):
                lines.append(str(params['gamma'][k][i][j]))

    # final three
    lines.extend([
        str(params['thickness_adjust']),
        str(params['metric_adjust']),
        str(int(params['restart']))
    ])

    in_file = os.path.join(output_dir, f'run_input_{sim_name}.txt')
    with open(in_file,'w') as f:
        f.write("\n".join(lines))
    return in_file


def plot_uv_surface(x_func, y_func, z_func, umin, umax, vmin, vmax, nu=50, nv=50):
    u = np.linspace(umin, umax, nu)
    v = np.linspace(vmin, vmax, nv)
    U, V = np.meshgrid(u, v)
    X = x_func(U, V)
    Y = y_func(U, V)
    Z = z_func(U, V)

    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    # compute the span of each axis
    # rx = X.max() - X.min()
    # ry = Y.max() - Y.min()
    # rz = Z.max() - Z.min()

    # set the box aspect to the actual data-ratios
    # ax.set_box_aspect((rx, ry, rz))
    # ax.set_box_aspect((1, 1, 1))

    
    surf = ax.plot_surface(X, Y, Z, 
                           rstride=1, cstride=1, edgecolor='none')
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    fig.colorbar(surf, ax=ax, shrink=0.5, label='Z value')
    plt.tight_layout()
    return ax


def func_to_matlab_str(func):
    """
    Turn a lambda like `lambda u,v: np.sin(u)*np.cos(v)` or `lambda u,v: u**2`
    into a MATLAB/Octave‐style string: 'sin(u)*cos(v)' or 'u^2', etc.
    """
    # grab the RHS of the lambda
    src = inspect.getsource(func).strip()
    m = re.match(r'.*lambda\s*[\w, ]*:\s*(.*)', src)
    if not m:
        raise ValueError("Expected a one-liner lambda(u,v): …")
    expr = m.group(1)

    # 3) strip off the "np." (or "numpy.") prefix everywhere
    expr = re.sub(r'\b(?:np|numpy)\.([A-Za-z_]\w*)', r'\1', expr)

    # 4) sympify with *all* Sympy names available
    sym = sp.sympify(expr, locals={**_sympy_namespace, 'u': u, 'v': v})

    # 5) print out Octave/MATLAB code
    return octave_code(sym)

def calc_Gamma(a, u, v):
    """
    give a 2x2 matrix a with sympy functions
    returns the Levi Chivita connectoin
    """
    inva  = a.inv()
    da_du = a.diff(u)
    da_dv = a.diff(v)

    Gamma = [[[sp.simplify(sp.Rational(1, 2) * sum(
                    inva[k, l] * (
                    (da_du[l, j] if i == 0 else da_dv[l, j]) +
                    (da_du[l, i] if j == 0 else da_dv[l, i]) -
                    (da_du[i, j] if l == 0 else da_dv[i, j])
                    )
                for l in range(2)))
            for j in range(2)]
            for i in range(2)]
            for k in range(2)]

    Gamma_str = [[[octave_code(Gamma[k][i][j])
                for j in range(2)]
                for i in range(2)]
                for k in range(2)]
    return Gamma_str