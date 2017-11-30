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
x is axial
"""
function beam_matrix(L, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz)

    # initialize
    K = zeros(2*DOF, 2*DOF)  # stiffness matrix
    M = zeros(2*DOF, 2*DOF)  # inertial matrix
    F = zeros(2*DOF)  # force vector

    # --- axial ----
    idx = [1, 7]

    kx = axialstiffness(EA, L)
    K[idx, idx] = kx

    mx = axialinertial(rhoA, L)
    M[idx, idx] = mx

    fx = axialloads(Px, L)
    F[idx] = fx

    # --- torsion ---
    idx = [2, 8]

    kt = axialstiffness(GJ, L)
    K[idx, idx] = kt

    mt = axialinertial(rhoJ, L)
    M[idx, idx] = mt

    # no distributed torsional loads

    # --- bending in y ---
    idx = [3, 4, 9, 10]

    ky = bendingstiffness(EIy, L)
    K[idx, idx] = ky

    my = bendinginertial(rhoA, L)
    M[idx, idx] = my

    fy = bendingloads(Py, L)
    F[idx] = fy


    # --- bending in z ---
    kz = bendingstiffness(EIy, L)
    idx = [5, 6, 11, 12]
    K[idx, idx] = kz

    mz = my  # inertial is same in y and z
    M[idx, idx] = mz

    fz = bendingloads(Pz, L)
    F[idx] = fz


    return K, M, F

end



# assembles FEA matrices for the various elements into global matrices for the structure
# EIx, EA etc. are arrays of length nodes
# matrix are of size DOF*nodes x DOF*nodes
# addedK: additional stiffness (or infinite stiffness for a rigid connection)
function fea_analysis(x, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
    Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    nodes = length(x)
    elements = nodes - 1

    # --- assemble global matrices -----
    K = zeros(DOF*nodes, DOF*nodes)
    M = zeros(DOF*nodes, DOF*nodes)
    F = zeros(DOF*nodes)

    for i = 1:elements
        start = (i-1)*DOF  # (0, 0) start of matrix

        Ksub, Msub, Fsub = beam_matrix(x[i+1] - x[i], EIy[i:i+1], EIz[i:i+1], EA[i:i+1],
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
        println("WARNING: structure improperly constrained if deflections are desired.")  #  (Exception: ", err, ")")
    end

    # ----- compute eigenvalues -----
    lambda, V = eig(K, M)  # eigenvalues are omega^2 (currently not using eigenvectors)

    # don't save rigid modes
    lambda = lambda[lambda .> 1e-6]

    # convert to freq
    freq = sqrt.(lambda) / (2*pi)

    return delta, freq
end


"""
Currently it assumes that the free end is at index 1.  TODO: allow either.
"""
function strain(x, y, z, EIy, EIz, EA, Px, Py, Pz, Fxpt, Fypt, Fzpt, Mxpt, Mypt, Mzpt)

    # initialize
    n = length(x)
    Nx = zeros(n)
    Vy = zeros(n)
    Vz = zeros(n)
    Tx = zeros(n)
    My = zeros(n)
    Mz = zeros(n)

    # integrate
    for i = 1:n-1
        Nx[i+1] = Nx[i] + Fxpt[i] + (x[i+1] - x[i])*(Px[i] + Px[i+1])/2.0
        Vy[i+1] = Vy[i] + Fypt[i] + (x[i+1] - x[i])*(Py[i] + Py[i+1])/2.0
        Vz[i+1] = Vz[i] + Fzpt[i] + (x[i+1] - x[i])*(Pz[i] + Pz[i+1])/2.0

        Tx[i+1] = Tx[i] + Mxpt[i]
        My[i+1] = My[i] + Mypt[i] - (x[i+1] - x[i])*(Vz[i] + Fzpt[i]) - (x[i+1] - x[i])^2*(2*Pz[i] + Pz[i+1])/6.0
        Mz[i+1] = Mz[i] + Mzpt[i] + (x[i+1] - x[i])*(Vy[i] + Fypt[i]) + (x[i+1] - x[i])^2*(2*Py[i] + Py[i+1])/6.0
    end

    # strain
    epsilon_axial = Tx./EA + Mz./EIz.*y - My./EIy.*z

    return epsilon_axial, Nx, Vy, Vz, Tx, My, Mz
end


end # module
