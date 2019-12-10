module BeamFEA

import LinearAlgebra: Symmetric, eigvals
import Interpolations: LinearInterpolation

export BeamProperties, Loads, NoLoads, Stiffness, Deflection
export fea, strain

const DOF = 6

"""
wrapper for linear interpolation
"""
function linterp(x, y, xpt)
    intp = LinearInterpolation(x, y)
    return intp(xpt)
end

"""
(private)
see doc: stiffness matrix
"""
function bendingstiffness(EI, L)

    e1 = EI[1]
    e2 = EI[2]

    k = [
        (6*e1 + 6*e2)/L^3   (4*e1 + 2*e2)/L^2   (-6*e1 - 6*e2)/L^3  (2*e1 + 4*e2)/L^2;
        0.0                 (3*e1 + e2)/L     (-4*e1 - 2*e2)/L^2  (e1 + e2)/L;
        0.0                 0.0                 (6*e1 + 6*e2)/L^3   (-2*e1 - 4*e2)/L^2;
        0.0                 0.0                 0.0                 (e1 + 3*e2)/L;
    ]

    return Symmetric(k)
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

    return Symmetric(k)
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
function bendinginertial(rhoA, L)

    a1 = rhoA[1]
    a2 = rhoA[2]

    m = 1.0/420*[
        (120*a1 + 36*a2)*L  (15*a1 + 7*a2)*L^2      (27*a1 + 27*a2)*L   (-7*a1 - 6*a2)*L^2;
        0.0                 (2.5*a1 + 1.5*a2)*L^3   (6*a1 + 7*a2)*L^2   (-1.5*a1 - 1.5*a2)*L^3;
        0.0                 0.0                     (36*a1 + 120*a2)*L  (-7*a1 - 15*a2)*L^2;
        0.0                 0.0                     0.0                 (1.5*a1 + 2.5*a2)*L^3;
    ]

    return Symmetric(m)
end


"""
(private)
see doc: Inertial Matrix
"""
function axialinertial(rhoA, L)

    r1 = rhoA[1]
    r2 = rhoA[2]

    m = L*[
        (1.0/4*r1 + 1.0/12*r2)   (1.0/12*r1 + 1.0/12*r2);
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
    type = typeof(EIy[1])
    K = zeros(type, 2*DOF, 2*DOF)  # stiffness matrix
    M = zeros(type, 2*DOF, 2*DOF)  # inertial matrix
    F = zeros(type, 2*DOF)  # force vector

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
    idx = [5, 6, 11, 12]
    
    kz = bendingstiffness(EIz, L)
    K[idx, idx] = kz

    M[idx, idx] = my  # inertial is same in y and z

    fz = bendingloads(Pz, L)
    F[idx] = fz


    return K, M, F
end


struct BeamProperties{TF, TA<:AbstractVector{TF}}

    x::TA  # sometimes this is a range
    EIy::Vector{TF}
    EIz::Vector{TF}
    EA::Vector{TF}
    GJ::Vector{TF}
    rhoA::Vector{TF}
    rhoJ::Vector{TF}

    # EIx, EA etc. are arrays of length nodes
end

struct Loads{TF, TA<:AbstractVector{TF}}

    x::TA
    Px::Vector{TF}
    Py::Vector{TF}
    Pz::Vector{TF}
    Fx::Vector{TF}
    Fy::Vector{TF}
    Fz::Vector{TF}
    Mx::Vector{TF}
    My::Vector{TF}
    Mz::Vector{TF}
end

NoLoads(x) = Loads(x, zeros(length(x)), zeros(length(x)), zeros(length(x)), zeros(length(x)), 
    zeros(length(x)), zeros(length(x)), zeros(length(x)), zeros(length(x)), zeros(length(x)))

# addedK: additional stiffness (or infinite stiffness for a rigid connection)
struct Stiffness{TF}
    
    kx::Vector{TF}
    ky::Vector{TF}
    kz::Vector{TF}
    kthetax::Vector{TF}
    kthetay::Vector{TF}
    kthetaz::Vector{TF}
end

struct Deflection{TF}
    x::Vector{TF}
    y::Vector{TF}
    z::Vector{TF}
    thetax::Vector{TF}
    thetay::Vector{TF}
    thetaz::Vector{TF}
end

struct ShearBending{TF}
    Nx::Vector{TF}
    Vy::Vector{TF}
    Vz::Vector{TF}
    Tx::Vector{TF}
    My::Vector{TF}
    Mz::Vector{TF}
end

"""
    fea_analysis(beam, loads, stiffness)

Beam finite element analysis.

**Inputs**
- beam::BeamProperties
- loads::Loads
- stiffness::Stiffness

**Outputs**
- defl::Deflection
- freq: natural frequencies
"""
function fea(beam, loads, stiffness)

    nodes = length(beam.x)
    elements = nodes - 1

    # ---- interpolate loads onto beam ----
    Px = linterp(loads.x, loads.Px, beam.x)
    Py = linterp(loads.x, loads.Py, beam.x)
    Pz = linterp(loads.x, loads.Pz, beam.x)
    Fx = linterp(loads.x, loads.Fx, beam.x)
    Fy = linterp(loads.x, loads.Fy, beam.x)
    Fz = linterp(loads.x, loads.Fz, beam.x)
    Mx = linterp(loads.x, loads.Mx, beam.x)
    My = linterp(loads.x, loads.My, beam.x)
    Mz = linterp(loads.x, loads.Mz, beam.x)
    # ----------------------------------

    # ---- global to fea convention ----
    My, Mz = -Mz, My
    # EIy/z and kthetay/z changes done in place
    # ---------------------------------

    # --- assemble global matrices -----
    type = typeof(beam.EIy[1])
    K = zeros(type, DOF*nodes, DOF*nodes)
    M = zeros(type, DOF*nodes, DOF*nodes)
    F = zeros(type, DOF*nodes)

    # pull out submatrix for each element
    for i = 1:elements
        
        Ksub, Msub, Fsub = beam_matrix(beam.x[i+1] - beam.x[i], 
            beam.EIz[i:i+1], beam.EIy[i:i+1], # swapped y and z (global -> fea convention)
            beam.EA[i:i+1], beam.GJ[i:i+1], beam.rhoA[i:i+1], beam.rhoJ[i:i+1], 
            Px[i:i+1], Py[i:i+1], Pz[i:i+1])
        
        start = (i-1)*DOF
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

        kvec = [stiffness.kx[i], stiffness.kthetax[i], stiffness.ky[i], stiffness.kthetaz[i],  # swapped y and z (global -> fea convention)
            stiffness.kz[i], stiffness.kthetay[i]]  # order for convenience

        for j = 1:DOF  # iterate through each DOF

            if kvec[j] == Inf  # if rigid, remove this index
                save[start+j] = false
            else  # otherwise add directly to stiffness matrix
                K[start+j, start+j] += kvec[j]
            end
        end
    end

    K = K[save, save]
    M = M[save, save]
    F = F[save]

    # denote symmeric more efficient linear solve and eigenanalysis
    K = Symmetric(K)
    M = Symmetric(M)

    # ---- compute deflections -----
    # initialize to zero so we can keep rigid deflections
    delta = zeros(DOF*nodes)

    try
        deltasub = K\F
        delta[save] = deltasub  # insert nonzero deflections
    catch err
        @warn "structure improperly constrained if deflections are desired."
    end

    defl = Deflection(delta[1:DOF:end], delta[3:DOF:end], delta[5:DOF:end],
        delta[2:DOF:end], delta[6:DOF:end], -delta[4:DOF:end])  # swapped y and z (fea -> global convention)

    # ----- compute eigenvalues -----
    lambda = eigvals(K, M)  # eigenvalues are omega^2 (currently not using eigenvectors)

    # don't save rigid modes
    lambda = lambda[lambda .> 1e-6]

    # convert to freq
    freq = sqrt.(lambda) / (2*pi)

    return defl, freq
end


"""
Currently it assumes that the free end is at index 1.  TODO: allow either.
"""
function strain(y, z, beam, loads)

    # initialize
    n = length(beam.x)
    type = typeof(beam.EIy[1])
    Nx = Vector{type}(undef, n)
    Vy = Vector{type}(undef, n)
    Vz = Vector{type}(undef, n)
    Tx = Vector{type}(undef, n)
    My = Vector{type}(undef, n)
    Mz = Vector{type}(undef, n)

    # unpack
    x = beam.x
    Px = loads.Px
    Py = loads.Py
    Pz = loads.Pz
    Fxpt = loads.Fx
    Fypt = loads.Fy
    Fzpt = loads.Fz
    Mxpt = loads.Mx
    Mypt = loads.My
    Mzpt = loads.Mz

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
    epsilon_axial = Tx./beam.EA - Mz./beam.EIz.*y + My./beam.EIy.*z

    return epsilon_axial, ShearBending(Nx, Vy, Vz, Tx, My, Mz)
end


end # module
