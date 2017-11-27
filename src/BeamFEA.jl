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
    (4*e1 + 2*e2)/L^2   (3*e1 + 1*e2)/L^1   (-4*e1 - 2*e2)/L^2   (e1 + e2)/L;
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
computes FEM matrices for one 12-dof beam element

q = [x1, thetax1, y1, thetay1, z1, thetaz1, ... repeat for 2]
z is axial
"""
function beam_matrix(L, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz)

    # ---- fill matrices ----
    K = zeros(2*DOF, 2*DOF)  # stiffness matrix
    F = zeros(2*DOF)  # force vector

    # axial
    kz = axialstiffness(EA, L)
    idx = [5, 11]
    K[idx, idx] = kz

    fz = axialloads(Pz, L)
    F[idx] = fz

    # torsion
    kt = axialstiffness(GJ, L)
    idx = [6, 12]
    K[idx, idx] = kt

    # bending in x
    kx = bendingstiffness(EIx, L)
    idx = [1, 2, 7, 8]
    K[idx, idx] = kx

    fx = bendingloads(Px, L)
    F[idx] = fx

    # bending in y
    ky = bendingstiffness(EIy, L)
    idx = [3, 4, 9, 10]
    K[idx, idx] = ky

    fy = bendingloads(Py, L)
    F[idx] = fy

    return K, F

end



# assembles FEA matrices for the various elements into global matrices for the structure
# EIx, EA etc. are arrays of length nodes
# matrix are of size DOF*nodes x DOF*nodes
# addedK: additional stiffness (or infinite stiffness for a rigid connection)
function fea_analysis(z, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
    kx, ky, kz, kthetax, kthetay, kthetaz)

    nodes = length(z)
    elements = nodes - 1

    # --- assemble global matrices -----
    K = zeros(DOF*nodes, DOF*nodes)
    F = zeros(DOF*nodes)

    for i = 1:elements

        Ksub, Fsub = beam_matrix(z[i+1] - z[i], EIx[i:i+1], EIy[i:i+1], EA[i:i+1],
            GJ[i:i+1], rhoA[i:i+1], rhoJ[i:i+1], Px[i:i+1], Py[i:i+1], Pz[i:i+1])

        start = (i-1)*2*DOF+1
        finish = i*2*DOF
        idx = start:finish
        K[idx, idx] = Ksub
        F[idx] = Fsub
    end

    # --- apply boundary conditions ----

    # add additional stiffness (and save nodes to remove)
    save = trues(DOF*nodes)  # nodes to keep
    for i = 1:nodes

        stiffness = [kx[i], kthetax[i], ky[i], kthetay[i], kz[i], kthetaz[i]]  # order for convenience

        for j = 1:DOF  # iterate through each DOF

            idx = (i-1)*DOF+j  # corresponding index
            if stiffness[i] == Inf  # if rigid, remove this index
                save[idx] = false
            else  # otherwise add directly to stiffness matrix
                K[idx, idx] += stiffness[i]
            end
        end
    end

    K = K[save, save]
    F = F[save]


    # ---- compute deflections
    deltasub = K\F

    # add back in rigid deflections (0)
    delta = zeros(DOF*nodes)
    delta[save] = deltasub

    return delta
end


end # module
