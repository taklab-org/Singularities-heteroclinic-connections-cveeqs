# Codes of "Singularities and heteroclinic connections in complex-valued evolutionary equations with a quadratic nonlinearity"

This repository contains the MATLAB codes associated with the paper:
"Singularities and heteroclinic connections in complex-valued evolutionary equations with a quadratic nonlinearity"
by J Jaquette, J-P Lessard and A Takayasu.

**Abstract**: In this paper, we consider the dynamics of solutions to complex-valued evolutionary partial differential equations (PDEs) and show existence of heteroclinic orbits from nontrivial equilibria to zero via computer-assisted proofs. We also show that the existence of unbounded solutions along unstable manifolds at the equilibrium follows from the existence of heteroclinic orbits. Our computer-assisted proof consists of three separate techniques of rigorous numerics: an enclosure of a local unstable manifold at the equilibria, a rigorous integration of PDEs, and a constructive validation of a trapping region around the zero equilibrium.

These codes require *MATLAB* with [*INTLAB* - INTerval LABoratory](http://www.ti3.tu-harburg.de/rump/intlab/) (MATLAB toolbox for interval arithmetic) version 11 and [*Chebfun* - numerical computing with functions](https://www.chebfun.org/) version 5.7.0.

---

A rough correspondence for some of the files & our computer-assisted proofs in the present paper are as follows:

### Proof of an eigenpar for the parameterization method

```
>> cd Proofs_Eigenpairs
>> script_verify_eigenpairs
```

For the choice of eigenvector one can see the specific rescaling in `F_eigenpairs.m`.

### Proof of heteroclinic orbits coming out of the 1 complex dimensional local unstable manifold

For the case of $\theta=0$

```
>> cd ../time_integration/
>> script_proof_heteroclinics_CGL_0
```

For the case of $\theta=\pi/4$

```
>> script_proof_heteroclinics_CGL_pi_4
```

For the case of $\theta=\pi/2$ (nonlinear Schrödinger case)

```
>> script_proof_heteroclinics_NLS
```

Note that the code above gives computer-assisted proofs of 360 heteroclinic orbits. Each proof takes _almost an hour_ of computation time. One can access the results of computer-assisted proofs in the folder `CGL_0_pt1/`.

```
>> cd ../CGL_0_pt1/
>> script_plot_sol_dist_from_data
```

This code plots Figure 5(b). Also one can draw pictures both Figure 6(b) and 7(b) by using

```
>> cd ../CGL_pi_4_pt1/
>> script_plot_sol_dist_from_data
>> cd ../NLS/
>> script_plot_sol_dist_from_data
```

---

Copyright (C) 2021  J Jaquette, J-P Lessard and A Takayasu.
