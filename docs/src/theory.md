# Theory

This theory mostly follows that described by Yang[^1], although we use linearly varying elements rather than constant elements.  

## Stiffness Matrix

We assume that stiffness properties vary linearly between nodes:

```math 
EI(\eta) = EI_1 + (EI_2 - EI_1) \eta
```

The stiffness matrix is given by integrating (5.13, Yang):
```math
\begin{aligned}
k_{ij} &= \int_0^L EI f_i^{\prime\prime} f_j^{\prime\prime} dx \\
  &= L\int_0^1 EI(\eta) f_i^{\prime\prime}(\eta) f_j^{\prime\prime}(\eta) d\eta \\
\end{aligned}
```

The second derivative of the shape functions for beams are given as (5.10, Yang):
```math
\begin{aligned}
f_1^{\prime\prime}(\eta) &= \frac{1}{L^2}\left( -6 + 12 \eta \right)\\
f_2^{\prime\prime}(\eta) &= \frac{1}{L} \left( -4 + 6 \eta\right)\\
f_3^{\prime\prime}(\eta) &= \frac{1}{L^2}\left( 6 - 12 \eta \right)\\
f_4^{\prime\prime}(\eta) &= \frac{1}{L} \left( -2 + 6 \eta\right)\\
\end{aligned}
```

The integrals are worked out below.

```julia
using Polynomials

f1pp = Poly([-6.0, 12.0])
f2pp = Poly([-4.0, 6.0])
f3pp = Poly([6.0, -12.0])
f4pp = Poly([-2.0, 6.0])

e1 = Poly([1.0, -1.0]) 
e2 = Poly([0.0, 1.0])

function coeff(f1, f2, L1, L2)

    pint = polyint(e1*f1*f2)
    c1 = polyval(pint, 1.0) - polyval(pint, 0.0)

    pint = polyint(e2*f1*f2)
    c2 = polyval(pint, 1.0) - polyval(pint, 0.0)

    println(c1, " ", c2, " ", L1 + L2 - 1)
    
end

L1 = 2.0
L2 = 1.0
L3 = 2.0
L4 = 1.0

coeff(f1pp, f1pp, L1, L1)
coeff(f1pp, f2pp, L1, L2)
coeff(f1pp, f3pp, L1, L3)
coeff(f1pp, f4pp, L1, L4)
println()

coeff(f2pp, f2pp, L2, L2)
coeff(f2pp, f3pp, L2, L3)
coeff(f2pp, f4pp, L2, L4)
println()

coeff(f3pp, f3pp, L3, L3)
coeff(f3pp, f4pp, L3, L4)
println()

coeff(f4pp, f4pp, L4, L4)
```

Which produces the following output:
```
6.0 6.0 3
4.0 2.0 2
-6.0 -6.0 3
2.0 4.0 2

3.0 1.0 1
-4.0 -2.0 2
1.0 1.0 1

6.0 6.0 3
-2.0 -4.0 2

1.0 3.0 1
```

Mathematically, we can express the (symmetric) stiffness matrix as:

```math
\begin{bmatrix}
    F_{y1} \\ M_{y1} \\ F_{y2} \\ M_{y2}
\end{bmatrix}
= 
\begin{bmatrix}
    \frac{6EI_1 + 6 EI_2}{L^3} & \frac{4EI_1 + 2 EI_2}{L^2} & \frac{-6EI_1 - 6 EI_2}{L^3} &\frac{2EI_1 + 4 EI_2}{L^2}\\
    sym & \frac{3EI_1 + 1 EI_2}{L} & \frac{-4EI_1 - 2 EI_2}{L^2} &\frac{EI_1 + EI_2}{L}\\
    sym & sym & \frac{6EI_1 + 6 EI_2}{L^3} &\frac{-2 EI_1  -4 EI_2}{L^2}\\
    sym & sym & sym &\frac{ EI_1 + 3 EI_2}{L}\\
\end{bmatrix}
\begin{bmatrix}
    y_1 \\ \theta_{y1} \\ y_2 \\ \theta_{y2}
\end{bmatrix}
```
and similarly in z.

For axial (and torsional) stiffness, the stiffness matrix is given by (4.9, Yang):
``k_{ij} = L \int_0^1 EA(\eta) f_i^\prime(\eta) f_j^\prime(eta) d\eta``
and the derivatives of the axial shape functions are:
```math
\begin{aligned}
f_1^\prime &= -1/L\\
f_2^\prime &= 1/L
\end{aligned}
```

We assume stiffness properties vary linearly yielding:

```math
\begin{bmatrix}
    F_{x1} \\ F_{x2}
\end{bmatrix}
= \frac{EA_1 + EA_2}{2 L}
\begin{bmatrix}
    1 & -1 \\
    -1 & 1
\end{bmatrix}
\begin{bmatrix}
    x_1 \\ x_2
\end{bmatrix}
```
Torsion is the same form:
```math
\begin{bmatrix}
    M_{x1} \\ M_{x2}
\end{bmatrix}
= \frac{GJ_1 + GJ_2}{2 L}
\begin{bmatrix}
    1 & -1 \\
    -1 & 1
\end{bmatrix}
\begin{bmatrix}
    \theta_{x1} \\ \theta_{x2}
\end{bmatrix}
```

## Distributed Loads

We assume that distributed loads vary linearly between nodes:
``p(\eta) = p_1 + (p_2 - p_1) \eta``

The equivalent concentrated loads are given by integrating (5.73, Yang):
```math
\begin{aligned}
F_i &= \int_0^L p(x) f_i(x) dx \\
&= L \int_0^1 p(\eta) f_i(\eta) d\eta \\
\end{aligned}
```

The shape functions for beams are given as (5.7, Yang):
```math
\begin{aligned}
f_1(\eta) &= 1 - 3\eta^2 + 2 \eta^3\\
f_2(\eta) &= L \left( \eta - 2\eta^2 + \eta^3 \right)\\
f_3(\eta) &= 3\eta^2 - 2 \eta^3\\
f_4(\eta) &= L \left( - \eta^2 + \eta^3 \right)
\end{aligned}
```
The integrals are worked out below.

```julia


f1 = Poly([1.0, 0.0, -3.0, 2.0])
f2 = Poly([0.0, 1.0, -2.0, 1.0])
f3 = Poly([0.0, 0.0, 3.0, -2.0])
f4 = Poly([0.0, 0.0, -1.0, 1.0])

p1 = Poly([1.0, -1.0])  # p2 is 0
p2 = Poly([0.0, 1.0])  # p1 is 0

function coeff(f)

    pint = polyint(p1*f)
    c1 = polyval(pint, 1.0) - polyval(pint, 0.0)

    pint = polyint(p2*f)
    c2 = polyval(pint, 1.0) - polyval(pint, 0.0)

    println(c1, " ", c2)
    
end

coeff(f1)
coeff(f2)
coeff(f3)
coeff(f4)
```

```
0.35 0.15000000000000002
0.050000000000000044 0.033333333333333326
0.15000000000000002 0.35
-0.033333333333333326 -0.04999999999999999
```

mathematically:
```math
\begin{aligned}
F_{y1} &= L \left[ \frac{7}{20} p_1 + \frac{3}{20} p_2 \right]\\
M_{y1} &= L^2 \left[ \frac{1}{20} p_1 + \frac{1}{30} p_2 \right]\\
F_{y2} &= L \left[ \frac{3}{20} p_1 + \frac{7}{20} p_2 \right]\\
M_{y2} &= L^2 \left[ -\frac{1}{30} p_1 - \frac{1}{20} p_2 \right]
\end{aligned}
```
and similarly for z.

The axial loads are similar, but simpler to calculate.  Using the shape functions for the axial loads:
```math
\begin{aligned}
f_1 &= 1 - \eta\\
f_2 &= \eta\\
\end{aligned}
```
The integral yields:
```math
\begin{aligned}
F_{x1} &= L \left[ \frac{p_1}{3} + \frac{p_2}{6} \right]\\
F_{x2} &= L \left[ \frac{p_1}{6} + \frac{p_2}{3} \right]\\
\end{aligned}
```
It is assumed that the distributed loads ``P_x`` are aligned with the center of twist (or close enough) so there are no distributed torsional loads.


## Inertial Matrix

We assume that mass properties vary linearly between nodes:
``\rho A(\eta) = \rho A_1 + (\rho A_2 - \rho A_1) \eta``

The mass matrix integral is the same for bending and the axial direction (7.18 and 7.29, Yang):
```math
\begin{aligned}
m_{ij} &= \int_0^L \rho A f_i f_j dx \\
  &= L\int_0^1 \rho A(\eta) f_i(\eta) f_j(\eta) d\eta \\
\end{aligned}
```

The shape functions were given previously under Distributed Loads and the integrals are worked out below.

```julia

f1 = Poly([1.0, 0.0, -3.0, 2.0])
f2 = Poly([0.0, 1.0, -2.0, 1.0])
f3 = Poly([0.0, 0.0, 3.0, -2.0])
f4 = Poly([0.0, 0.0, -1.0, 1.0])

A1 = Poly([1.0, -1.0])  # A2 is 0
A2 = Poly([0.0, 1.0])  # A1 is 0

function coeff(f1, f2, L1, L2)

    pint = polyint(A1*f1*f2)
    c1 = polyval(pint, 1.0) - polyval(pint, 0.0)

    pint = polyint(A2*f1*f2)
    c2 = polyval(pint, 1.0) - polyval(pint, 0.0)

    println(c1*420, " ", c2*420, " ", L1 + L2 + 1)
    
end

L1 = 0
L2 = 1
L3 = 0
L4 = 1

coeff(f1, f1, L1, L1)
coeff(f1, f2, L1, L2)
coeff(f1, f3, L1, L3)
coeff(f1, f4, L1, L4)
println()

coeff(f2, f2, L2, L2)
coeff(f2, f3, L2, L3)
coeff(f2, f4, L2, L4)
println()

coeff(f3, f3, L3, L3)
coeff(f3, f4, L3, L4)
println()

coeff(f4, f4, L4, L4)
```

```
119.99999999999994 36.00000000000006 1
15.000000000000039 6.999999999999952 2
27.000000000000092 26.99999999999995 1
-6.999999999999952 -6.000000000000002 2

2.499999999999921 1.4999999999999947 3
5.999999999999885 7.000000000000021 2
-1.4999999999999947 -1.4999999999999947 3

35.999999999999964 120.00000000000004 1
-7.000000000000021 -14.999999999999993 2

1.4999999999999947 2.5000000000000027 3
```

Mathematically, we can express the (symmetric) inertial matrix as:

```math
\begin{bmatrix}
    F_{y1} \\ M_{y1} \\ F_{y2} \\ M_{y2}
\end{bmatrix}
= 
\frac{1}{420}
\begin{bmatrix}
    (120 \rho A_1 + 36 \rho A_2) L & (15 \rho A_1 + 7 \rho A_2) L^2 & (27 \rho A_1 + 27 \rho A_2) L & (-7 \rho A_1 -6  \rho A_2) L^2\\
    sym & (2.5 \rho A_1 + 1.5 \rho A_2) L^3 & (6 \rho A_1 + 7 \rho A_2) L^2 & (-1.5 \rho A_1 - 1.5 \rho A_2) L^3\\
    sym & sym & (36 \rho A_1 + 120 \rho A_2) L & (-7 \rho A_1 -15 \rho A_2) L^2\\
    sym & sym & sym & (1.5 \rho A_1 + 2.5 \rho A_2) L^3\\
\end{bmatrix}
\begin{bmatrix}
    \ddot{y_1} \\ \ddot{\theta_{y1}} \\ y_2 \\ \theta_{y2}
\end{bmatrix}
```
and similarly for z.

In the axial direction the integral is the same but the shape functions are simpler:
```math
\begin{aligned}
f_1 &= 1 - \eta\\
f_2 &= \eta\\
\end{aligned}
```

```math
\begin{aligned}
m_{ij} = L\int_0^1 \rho A(\eta) f_i(\eta) f_j(\eta) d\eta \\
\end{aligned}
```

```julia
f1 = Poly([1.0, -1.0])
f2 = Poly([0.0, 1.0])

A1 = Poly([1.0, -1.0])  # A2 is 0
A2 = Poly([0.0, 1.0])  # A1 is 0

function coeff(f1, f2)

    pint = polyint(A1*f1*f2)
    c1 = polyval(pint, 1.0) - polyval(pint, 0.0)

    pint = polyint(A2*f1*f2)
    c2 = polyval(pint, 1.0) - polyval(pint, 0.0)

    println(c1, " ", c2)
    
end


coeff(f1, f1)
coeff(f1, f2)
coeff(f2, f2)
```

```
0.25 0.08333333333333337
0.08333333333333337 0.08333333333333331
0.08333333333333331 0.25
```

The result is:
```math
\begin{bmatrix}
    F_{x1} \\ F_{x2}
\end{bmatrix}
= L
\begin{bmatrix}
    \frac{1}{4}\rho A_1 + \frac{1}{12}\rho A_2 & \frac{1}{12}\rho A_1 + \frac{1}{12}\rho A_2 \\
    \frac{1}{12}\rho A_1 + \frac{1}{12}\rho A_2 & \frac{1}{12}\rho A_1 + \frac{1}{4}\rho A_2
\end{bmatrix}
\begin{bmatrix}
    \ddot{x_1} \\ \ddot{x_2}
\end{bmatrix}
```

Torsional inertia is of the same form:
The result is:
```math
\begin{bmatrix}
    M_{x1} \\ M_{x2}
\end{bmatrix}
= L
\begin{bmatrix}
    \frac{1}{4}\rho J_1 + \frac{1}{12}\rho J_2 & \frac{1}{12}\rho J_1 + \frac{1}{12}\rho J_2 \\
    \frac{1}{12}\rho J_1 + \frac{1}{12}\rho J_2 & \frac{1}{12}\rho J_1 + \frac{1}{4}\rho J_2
\end{bmatrix}
\begin{bmatrix}
    \ddot{\theta_{x1}} \\ \ddot{\theta_{x2}}
\end{bmatrix}
```

## Displacement

Displacement is easily computed by solving the linear system (after applying boundary conditions): 
``K u = F``

## Natural frequencies

Natural frequencies are found by solving the generalized eigenvalue problem (7.24, Yang)
``[K - \omega^2 M]u = 0``
Solving the generalized eigenvalue problem gives ``\lambda = \omega^2``.  We find the corresponding frequencies as:
``f = \frac{\sqrt{\lambda}}{2\pi}``

## Axial Stress

The computation of axial stress is actually independent from the finite element analysis and depends only on the loading  and geometry.  First the forces and moments must be integrated along the beam starting from a free-end where the forces and moments are known (they are zero unless there is an external point load right at the tip).  Each segment of the beam must be integrated sequentially.  

A given beam segment is defined below:

<img src="element.png" width="400px"/>

We assume that the distributed load varies linear, but separate distributions can exist in each direction (including axial in ``x``).  External point forces and moments ``F_{pt}`` and ``M_{pt}`` can exist at each node in all three directions.  From this definition we can integrate to find the forces and moments throughout the beam.  Positive moments we define for moments on the left side of the beam using the right-hand rule according to the defined coordinate system.

```math
\begin{aligned}
    V_{i+1} &= V_i + {F_{pt}}_{i} + (x_{i+1}-x_i) \int_{0}^1 q(\eta) d\eta \\
    M_{i+1} &= M_i + {M_{pt}}_{i} + (x_{i+1}-x_i) \int_{0}^1 V(\eta) d\eta
\end{aligned}
```

If the distributed load is linear:
``q(\eta) = q_i + \eta(q_{i+1} - q_i) ``
then:
```math
\begin{aligned}
    V(\eta) &= V_i + {F_{pt}}_{i} + (x_{i+1}-x_i) \int \left(q_i + \eta(q_{i+1} - q_i)\right) d\eta \\
    &= V_i + {F_{pt}}_{i} + (x_{i+1}-x_i) \left(q_i \eta + \frac{\eta^2}{2}(q_{i+1} - q_i)\right)\\
    M(\eta) &= M_i + {M_{pt}}_{i} + (x_{i+1}-x_i) \int \left( V_i + {F_{pt}}_{i} + (x_{i+1}-x_i) \left(q_i \eta + \frac{\eta^2}{2}(q_{i+1} - q_i)\right) \right) d\eta\\
    &= M_i + {M_{pt}}_{i} + (x_{i+1}-x_i)(V_i + {F_{pt}}_{i}) \eta +    (x_{i+1}-x_i)^2 \left(q_i \frac{\eta^2}{2} + \frac{\eta^3}{6}(q_{i+1} - q_i)\right) 
\end{aligned}
```

Evaluating at the end point (``\eta = 1``) and simplifying yields:
```math
\begin{aligned}
    V_{i+1} &= V_i + {F_{pt}}_{i} + (x_{i+1}-x_i) \left(\frac{q_i + q_{i+1}}{2} \right) \\
    M_{i+1} &= M_i + {M_{pt}}_{i} + (x_{i+1}-x_i)(V_i + {F_{pt}}_{i}) +  (x_{i+1}-x_i)^2 \left(\frac{2 q_i + q_{i+1}}{6}\right)
\end{aligned}
```

This applies for all three forces: ``V_x, V_y, N_z``.  For the moments, because this is a beam and we assume any axial loads act along the axis, the torsional loads ``T_z`` will not change excepting point torsion loads.  However, the bending moments ``M_y, M_z`` are affected by the distributed loads ``q_y`` and ``q_z``.  Explicitly using our coordinate sign conventions:
```math
\begin{aligned}
    {N_x}_{i+1} &= {N_x}_i + {{F_x}_{pt}}_{i} + (x_{i+1}-x_i) \left(\frac{{q_x}_i + {q_x}_{i+1}}{2} \right) \\
    {V_y}_{i+1} &= {V_y}_i + {{F_y}_{pt}}_{i} + (x_{i+1}-x_i) \left(\frac{{q_y}_i + {q_y}_{i+1}}{2} \right) \\
    {V_z}_{i+1} &= {V_z}_i + {{F_z}_{pt}}_{i} + (x_{i+1}-x_i) \left(\frac{{q_z}_i + {q_z}_{i+1}}{2} \right) \\
    {T_x}_{i+1} &= {T_x}_i + {{T_x}_{pt}}_{i} \\
    {M_y}_{i+1} &= {M_y}_i + {{M_y}_{pt}}_{i} - (x_{i+1}-x_i)({V_z}_i + {{F_z}_{pt}}_{i}) - (x_{i+1}-x_i)^2 \left(\frac{2 {q_z}_i + {q_z}_{i+1}}{6}\right)\\
    {M_z}_{i+1} &= {M_z}_i + {{M_z}_{pt}}_{i} + (x_{i+1}-x_i)({V_y}_i + {{F_y}_{pt}}_{i}) + (x_{i+1}-x_i)^2 \left(\frac{2 {q_y}_i + {q_y}_{i+1}}{6}\right)
\end{aligned}
```
With known shear/moment distribution and stiffness properties we can compute the stress as follows (or use ``E(y, z) = 1`` to compute strain):

```math
\sigma_{xx}(y, z) = E(y, z) \left(\frac{M_z}{[EI]_z} y - \frac{M_y}{[EI]_y} z + \frac{N_x}{[EA]} \right)
```


[^1]: Yang, T. Y., Finite element structural analysis, Prentice Hall, 1986.