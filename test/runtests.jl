using BeamFEA
using Base.Test

@testset "Simple Beam Deflection Tests" begin

    # --- Test data from "Finite Element Structural Analysis", Yang, pg. 145 ----
    # cantilever_deflection

    # arbitrary values
    E = 2.0
    I = 3.0
    L = 4.0
    p0 = 5.0

    z = [0, L]
    EIx = [E*I, E*I]
    EIy = [E*I, E*I]
    EA = [E*1.0, E*1.0]
    GJ = [1.0, 1.0]
    rhoA = [1.0, 1.0]
    rhoJ = [1.0, 1.0]
    Px = [0.0, 0.0]
    Py = [0.0, -p0]
    Pz = [0.0, 0.0]
    Fx = [0.0, 0.0]
    Fy = [0.0, 0.0]
    Fz = [0.0, 0.0]
    Mx = [0.0, 0.0]
    My = [0.0, 0.0]
    Mz = [0.0, 0.0]
    kx = [0.0, Inf]
    ky = [0.0, Inf]
    kz = [0.0, Inf]
    kthetax = [0.0, Inf]
    kthetay = [0.0, Inf]
    kthetaz = [0.0, Inf]

    delta, freq = BeamFEA.fea_analysis(z, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
        Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    @test delta[3] ≈ -p0*L^4/(30*E*I)  atol=1e-6
    @test delta[4] ≈ p0*L^3/(24*E*I)  atol=1e-6
    @test delta[1] ≈ 0.0  atol=1e-6
    @test delta[2] ≈ 0.0  atol=1e-6
    @test delta[5] ≈ 0.0  atol=1e-6
    @test delta[6] ≈ 0.0  atol=1e-6


    # --- Test data from "Finite Element Structural Analysis", Yang, pg. 180 ----
    # tapered beam

    # arbitrary values
    E = 2.0
    I = 3.0
    L = 4.0
    P = 5.0

    z = [0, L/2, L]
    EIx = [E*I, 5*E*I, 9*E*I]
    EIy = [E*I, 5*E*I, 9*E*I]
    EA = [1.0, 1.0, 1.0]
    GJ = [1.0, 1.0, 1.0]
    rhoA = [1.0, 1.0, 1.0]
    rhoJ = [1.0, 1.0, 1.0]
    Px = [0.0, 0.0, 0.0]
    Py = [0.0, 0.0, 0.0]
    Pz = [0.0, 0.0, 0.0]
    Fx = [0.0, 0.0, 0.0]
    Fy = [-P, 0.0, 0.0]
    Fz = [0.0, 0.0, 0.0]
    Mx = [0.0, 0.0, 0.0]
    My = [0.0, 0.0, 0.0]
    Mz = [0.0, 0.0, 0.0]
    kx = [0.0, 0.0, Inf]
    ky = [0.0, 0.0, Inf]
    kz = [0.0, 0.0, Inf]
    kthetax = [0.0, 0.0, Inf]
    kthetay = [0.0, 0.0, Inf]
    kthetaz = [0.0, 0.0, Inf]

    delta, freq = BeamFEA.fea_analysis(z, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
        Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    # println(delta)
    # println(-0.051166*P*L^3/(E*I))
    # println(0.090668*P*L^2/(E*I))
    @test delta[3] ≈ -0.051166*P*L^3/(E*I)  atol=5e-3
    @test delta[4] ≈ 0.090668*P*L^2/(E*I)  atol=1e-2
    @test delta[1] ≈ 0.0  atol=1e-6
    @test delta[2] ≈ 0.0  atol=1e-6
    @test delta[5] ≈ 0.0  atol=1e-6
    @test delta[6] ≈ 0.0  atol=1e-6

end  # end unit test


@testset "Beam Frequency Tests" begin


    # Test data from "Consistent Mass Matrix for Distributed Mass Systmes", John Archer,
    # Journal of the Structural Division Proceedings of the American Society of Civil Engineers,
    # pg. 168
    E = 2.0
    I = 3.0
    L = 4.0
    A = 5.0
    rho = 6.0

    n = 1
    nodes = n + 1
    z = linspace(0, L, nodes)
    EIx = E*I*ones(nodes)
    EIy = E*I*ones(nodes)
    EA = ones(nodes)
    GJ = ones(nodes)
    rhoA = rho*A*ones(nodes)
    rhoJ = ones(nodes)

    # TODO: add helper functions to initialize these to zero when not used
    Px = zeros(nodes)
    Py = zeros(nodes)
    Pz = zeros(nodes)
    Fx = zeros(nodes)
    Fy = zeros(nodes)
    Fz = zeros(nodes)
    Mx = zeros(nodes)
    My = zeros(nodes)
    Mz = zeros(nodes)
    kx = zeros(nodes)
    ky = zeros(nodes)
    kz = zeros(nodes)
    kthetax = zeros(nodes)
    kthetay = zeros(nodes)
    kthetaz = zeros(nodes)

    delta, freq = BeamFEA.fea_analysis(z, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
        Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    alpha = rho * A * (n*L)^4 / (840.0 * E * I)

    @test isapprox(freq[2], sqrt(0.85714 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[3], sqrt(0.85714 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[5], sqrt(10.0 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[6], sqrt(10.0 / alpha) / (2*pi), atol=1e-6)

end
