module BeamFEA

"""
(private)
see doc: stiffness matrix
"""
function bendingstiffness(EI::Array{Float64}(2), L::Float64)

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
computes FEM matrices for one 12-dof beam element

q = [x1, thetax1, y1, thetay1, z1, thetaz1, ... repeat for 2]
z is axial
"""
function beam_matrix(L, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz)

    DOF = 12

    # ---- fill matrices ----
    K = zeros(DOF, DOF)  # stiffness matrix
    F = zeros(DOF)  # force vector

    # axial
    kz = axialstiffness(EA, L)
    idx = [5, 11]
    K[idx, idx] = kz

    fz = bendingloads(Pz, L)
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

end




end # module
