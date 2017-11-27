using BeamFEA
using Base.Test

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

delta = BeamFEA.fea_analysis(z, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
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

delta = BeamFEA.fea_analysis(z, EIx, EIy, EA, GJ, rhoA, rhoJ, Px, Py, Pz,
    Fx, Fy, Fz, Mx, My, Mz, kx, ky, kz, kthetax, kthetay, kthetaz)

println(delta)
# @test delta[3] ≈ -p0*L^4/(30*E*I)  atol=1e-6
# @test delta[4] ≈ p0*L^3/(24*E*I)  atol=1e-6
# @test delta[1] ≈ 0.0  atol=1e-6
# @test delta[2] ≈ 0.0  atol=1e-6
# @test delta[5] ≈ 0.0  atol=1e-6
# @test delta[6] ≈ 0.0  atol=1e-6



println("tests passed!")
