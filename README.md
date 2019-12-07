# BeamFEA

<!--
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://byuflowlab.github.io/BeamFEA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://byuflowlab.github.io/BeamFEA.jl/dev)
-->

*A 6DOF beam finite element model*

Author: Andrew Ning

This is a finite element code for beam-like structures. The methodology uses Euler-Bernoulli beam elements with 12 degrees of freedom (3 translation and 3 rotational at each end of the element). Structural properties can vary linearly across each element.  BeamFEA can estimate structural mass, deflections in all degrees of freedom, coupled natural frequencies, critical global axial buckling loads, and axial stress/strain. All inputs and outputs are given about the elastic center and in principal axes in order to remove cross-coupling terms. 
