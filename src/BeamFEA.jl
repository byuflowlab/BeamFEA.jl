module BeamFEA

const DOF = 6


"""
(private)
see doc: stiffness matrix
"""
function bendingstiffness(EI::Array{Float64, 1}, L::Float64)

    e1 = EI[1]
    e2 = EI[2]

    k = [
    (6*e1 + 6*e2)/L^3   (4*e1 + 2*e2)/L^2   (-6*e1 - 6*e2)/L^3   (2*e1 + 4*e2)/L^2;
    (4*e1 + 2*e2)/L^2   (3*e1 + 1*e2)/L     (-4*e1 - 2*e2)/L^2   (e1 + e2)/L;
    (-6*e1 - 6*e2)/L^3  (-4*e1 - 2*e2)/L^2  (6*e1 + 6*e2)/L^3    (-2*e1 - 4*e2)/L^2;
    (2*e1 + 4*e2)/L^2   (e1 + e2)/L         (-2*e1 - 4*e2)/L^2   (1*e1 + 3*e2)/L;
    ]

    return k
end


"""
(private)
see doc: stiffness matrix
"""
function axialstiffness(EA, L)


    k = (EA[1] + EA[2])/(2*L)*[
            1.0 -1.0;
            -1.0  1.0
        ]

    return k
end

"""
(private)
see doc: Distributed Loads
"""
function bendingloads(P, L)

    p1 = P[1]
    p2 = P[2]

    f = L*[
    7.0/20*p1 + 3.0/20*p2;
    L*(1.0/20*p1 + 1.0/30*p2);
    3.0/20*p1 + 7.0/20*p2;
    L*(-1.0/30*p1 - 1.0/20*p2);
    ]

    return f
end


"""
(private)
see doc: Distributed Loads
"""
function axialloads(P, L)

    p1 = P[1]
    p2 = P[2]

    f = L*[
    p1/3.0 + p2/6.0;
    p1/6.0 + p2/3.0;
    ]

    return f
end


"""
(private)
see doc: Inertial Matrix
"""
function bendinginertial(rhoA::Array{Float64, 1}, L::Float64)

    a1 = rhoA[1]
    a2 = rhoA[2]

    m = 1.0/420*[
    (120*a1 + 36*a2)*L  (15*a1 + 7*a2)*L^2      (27*a1 + 27*a2)*L    (-7*a1 - 6*a2)*L^2;
    (15*a1 + 7*a2)*L^2  (2.5*a1 + 1.5*a2)*L^3   (6*a1 + 7*a2)*L^2    (-1.5*a1 - 1.5*a2)*L^3;
    (27*a1 + 27*a2)*L   (6*a1 + 7*a2)*L^2       (36*a1 + 120*a2)*L   (-7*a1 - 15*a2)*L^2;
    (-7*a1 - 6*a2)*L^2  (-1.5*a1 - 1.5*a2)*L^3  (-7*a1 - 15*a2)*L^2  (1.5*a1 + 2.5*a2)*L^3;
    ]

    return m
end


"""
(private)
see doc: Inertial Matrix
"""
function axialinertial(rhoA, L)

    r1 = rhoA[1]
    r2 = rhoA[2]

    m = L*[
    (1.0/4*r1 + 1.0/12*r2)   (1.0/12*r1 + 1.0/12*r2)
    (1.0/12*r1 + 1.0/12*r2)  (1.0/12*r1 + 1.0/4*r2)
    ]

    return m
end


"""
(private)
computes FEM matrices for one 12-dof beam element

q = [x1, thetax1, y1, thetay1, z1, thetaz1, ... repeat for 2]
z is axial
"""
function beam_matrix(L, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz)

    # initialize
    K = zeros(2*DOF, 2*DOF)  # stiffness matrix
    M = zeros(2*DOF, 2*DOF)  # inertial matrix
    F = zeros(2*DOF)  # force vector

    # --- axial ----
    idx = [5, 11]

    kz = axialstiffness(EA, L)
    K[idx, idx] = kz

    mz = axialinertial(rhoA, L)
    M[idx, idx] = mz

    fz = axialloads(Pz, L)
    F[idx] = fz

    # --- torsion ---
    idx = [6, 12]

    kt = axialstiffness(GJ, L)
    K[idx, idx] = kt

    mt = axialinertial(rhoJ, L)
    M[idx, idx] = mt

    # no distributed torsional loads

    # --- bending in x ---
    idx = [1, 2, 7, 8]

    kx = bendingstiffness(EIx, L)
    K[idx, idx] = kx

    mx = bendinginertial(rhoA, L)
    M[idx, idx] = mx

    fx = bendingloads(Px, L)
    F[idx] = fx


    # --- bending in y ---
    ky = bendingstiffness(EIy, L)
    idx = [3, 4, 9, 10]
    K[idx, idx] = ky

    my = mx  # inertial is same in x and y
    M[idx, idx] = my

    fy = bendingloads(Py, L)
    F[idx] = fy


    return K, M, F

end



# assembles FEA matrices for the various elements into global matrices for the structure
# EIx, EA etc. are arrays of length nodes
# matrix are of size DOF*nodes x DOF*nodes
# addedK: additional stiffness (or infinite stiffness for a rigid connection)
function fea_analysis(z, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
    Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    nodes = length(z)
    elements = nodes - 1

    # --- assemble global matrices -----
    K = zeros(DOF*nodes, DOF*nodes)
    M = zeros(DOF*nodes, DOF*nodes)
    F = zeros(DOF*nodes)

    for i = 1:elements
        start = (i-1)*DOF  # (0, 0) start of matrix

        Ksub, Msub, Fsub = beam_matrix(z[i+1] - z[i], EIx[i:i+1], EIy[i:i+1], EA[i:i+1],
            GJ[i:i+1], rhoA[i:i+1], rhoJ[i:i+1], Px[i:i+1], Py[i:i+1], Pz[i:i+1])

        idx = start+1:start+2*DOF
        K[idx, idx] += Ksub
        M[idx, idx] += Msub
        F[idx] += Fsub
    end

    # ---- add point loads ----
    for i = 1:nodes
        start = (i-1)*DOF

        F[start + 1] += Fx[i]
        F[start + 2] += Mx[i]
        F[start + 3] += Fy[i]
        F[start + 4] += My[i]
        F[start + 5] += Fz[i]
        F[start + 6] += Mz[i]
    end

    # ---- add point inertial elements -----


    # --- apply boundary conditions ----

    # add additional stiffness (and save nodes to remove)
    save = trues(DOF*nodes)  # nodes to keep
    for i = 1:nodes

        start = (i-1)*DOF

        stiffness = [kx[i], kthetax[i], ky[i], kthetay[i], kz[i], kthetaz[i]]  # order for convenience

        for j = 1:DOF  # iterate through each DOF

            if stiffness[j] == Inf  # if rigid, remove this index
                save[start+j] = false
            else  # otherwise add directly to stiffness matrix
                K[start+j, start+j] += stiffness[j]
            end
        end
    end

    K = K[save, save]
    M = M[save, save]
    F = F[save]


    # ---- compute deflections -----
    # initialize to zero so we can keep rigid deflections
    delta = zeros(DOF*nodes)

    try
        deltasub = K\F
        delta[save] = deltasub  # insert nonzero deflections
    catch err
        println("WARNING: structure improperly constrained.  Deflections incorrect. (Exception: ", err, ")")
    end

    # ----- compute eigenvalues -----
    lambda, V = eig(K, M)  # eigenvalues are omega^2 (currently not using eigenvectors)

    # don't save rigid modes
    lambda = lambda[lambda .> 1e-6]

    # convert to freq
    freq = sqrt.(lambda) / (2*pi)

    return delta, freq
end


end # module
