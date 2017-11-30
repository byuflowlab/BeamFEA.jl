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

    x = [0, L]
    EIy = [E*I, E*I]
    EIz = [E*I, E*I]
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

    delta, freq = BeamFEA.fea_analysis(x, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
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

    x = [0, L/2, L]
    EIy = [E*I, 5*E*I, 9*E*I]
    EIz = [E*I, 5*E*I, 9*E*I]
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

    delta, freq = BeamFEA.fea_analysis(x, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
        Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

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
    x = linspace(0, n*L, nodes)
    EIy = E*I*ones(nodes)
    EIz = E*I*ones(nodes)
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

    delta, freq = BeamFEA.fea_analysis(x, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
        Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    alpha = rho * A * (n*L)^4 / (840.0 * E * I)

    @test isapprox(freq[2], sqrt(0.85714 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[3], sqrt(0.85714 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[5], sqrt(10.0 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[6], sqrt(10.0 / alpha) / (2*pi), atol=1e-6)

    # simply supported
    kx[1] = Inf
    ky[1] = Inf
    kz[1] = Inf
    ky[end] = Inf
    delta, freq = BeamFEA.fea_analysis(x, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
        Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    @test isapprox(freq[2], sqrt(0.14286 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[5], sqrt(3.0 / alpha) / (2*pi), atol=1e-6)

    # second test from same paper (n = 2)

    n = 2
    nodes = n + 1
    x = linspace(0, n*L, nodes)
    EIy = E*I*ones(nodes)
    EIz = E*I*ones(nodes)
    EA = ones(nodes)
    GJ = ones(nodes)
    rhoA = rho*A*ones(nodes)
    rhoJ = ones(nodes)

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

    delta, freq = BeamFEA.fea_analysis(x, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
        Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    alpha = rho * A * (n*L)^4 / (840.0 * E * I)

    @test isapprox(freq[2], sqrt(0.59858 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[3], sqrt(0.59858 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[6], sqrt(5.8629 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[7], sqrt(5.8629 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[9], sqrt(36.659 / alpha) / (2*pi), atol=2e-6)
    @test isapprox(freq[10], sqrt(36.659 / alpha) / (2*pi), atol=2e-6)
    @test isapprox(freq[11], sqrt(93.566 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[12], sqrt(93.566 / alpha) / (2*pi), atol=1e-6)

    # simply supported
    kx[1] = Inf
    ky[1] = Inf
    kz[1] = Inf
    ky[end] = Inf
    delta, freq = BeamFEA.fea_analysis(x, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
        Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    @test isapprox(freq[2], sqrt(0.11688 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[5], sqrt(2.2858 / alpha) / (2*pi), atol=1e-6)
    @test isapprox(freq[8], sqrt(14.441 / alpha) / (2*pi), atol=1e-5)
    @test isapprox(freq[11], sqrt(48 / alpha) / (2*pi), atol=1e-6)


    # third test from same paper (n = 3)

    n = 3
    nodes = n + 1
    x = linspace(0, n*L, nodes)
    EIy = E*I*ones(nodes)
    EIz = E*I*ones(nodes)
    EA = ones(nodes)
    GJ = ones(nodes)
    rhoA = rho*A*ones(nodes)
    rhoJ = ones(nodes)

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

    delta, freq = BeamFEA.fea_analysis(x, EIy, EIz, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
        Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

    alpha = rho * A * (n*L)^4 / (840.0 * E * I)


    @test isapprox(freq[2], sqrt(0.59919 / alpha)/(2*pi), atol=1e-6)
    @test isapprox(freq[3], sqrt(0.59919 / alpha)/(2*pi), atol=1e-6)

    @test isapprox(freq[6], sqrt(4.5750 / alpha)/(2*pi), atol=1e-6)
    @test isapprox(freq[7], sqrt(4.5750 / alpha)/(2*pi), atol=1e-6)

    @test isapprox(freq[9], sqrt(22.010 / alpha)/(2*pi), atol=1e-6)
    @test isapprox(freq[10], sqrt(22.010 / alpha)/(2*pi), atol=1e-6)

    @test isapprox(freq[12], sqrt(70.920 / alpha)/(2*pi), atol=1e-6)
    @test isapprox(freq[13], sqrt(70.920 / alpha)/(2*pi), atol=1e-6)

    @test isapprox(freq[15], sqrt(265.91 / alpha)/(2*pi), atol=1e-5)
    @test isapprox(freq[16], sqrt(265.91 / alpha)/(2*pi), atol=1e-5)

    @test isapprox(freq[17], sqrt(402.40 / alpha)/(2*pi), atol=1e-6)
    @test isapprox(freq[18], sqrt(402.40 / alpha)/(2*pi), atol=1e-6)



end



@testset "shear and bending" begin

    # Test data from "Mechanical of Materials", Gere, 6th ed., pg. 273
    # cantilevered beam with linear distributed load

    L = 10.0
    q0 = 3.0

    n = 1
    nodes = n+1

    x = linspace(0, L, nodes)
    y = zeros(nodes)
    z = zeros(nodes)
    EIy = ones(nodes)
    EIz = ones(nodes)
    EA = ones(nodes)
    Px = [0.0, 0.0]
    Py = [0.0, -q0]
    Pz = [0.0, 0.0]
    Fx = zeros(nodes)
    Fy = zeros(nodes)
    Fz = zeros(nodes)
    Mx = zeros(nodes)
    My = zeros(nodes)
    Mz = zeros(nodes)

    epsilon, Nx, Vy, Vz, Tx, My, Mz = BeamFEA.strain(x, y, z, EIy, EIz, EA, Px, Py, Pz, Fx, Fy, Fz, Mx, My, Mz)

    @test isapprox(Vy[1], 0.0, atol=1e-6)
    @test isapprox(Vy[2], -q0*L/2.0, atol=1e-6)
    @test isapprox(Mz[1], 0.0, atol=1e-6)
    @test isapprox(Mz[2], -q0*L^2/6, atol=1e-6)


end
