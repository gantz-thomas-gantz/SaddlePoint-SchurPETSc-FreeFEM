# **Saddle-Point Systems and Schur Complement Preconditioning with PETSc and FreeFEM**

- **Authors** : Jan ESQUIVEL MARXEN, Thomas GANTZ

## Overview

This project addresses the numerical solution of **coupled partial differential equations (PDEs)** using the **finite element method**, in particular:

**A.** The **Poisson equation in 1D** with homogeneous Dirichlet boundary conditions and a volume constraint

> Find $u \in H^1_0(0,1)$ such that
>
> $$
> \begin{cases}
> -u''(x) = f(x) & \text{in } (0,1), \\
> \int_0^1 u(x)\,dx = C. & \text{(volume constraint)}
> \end{cases}
> $$

**B.** The **Stokes equations in 2D**

> Find $\boldsymbol{u} \in [H^1(\Omega)]^2$, $p \in L^2(\Omega)$ such that
>
> $$
> \begin{cases}
> -\nu \Delta \boldsymbol{u} + \nabla p = \boldsymbol{f} & \text{in } \Omega, \\
> \text{div}(\boldsymbol{u}) = 0 & \text{in } \Omega, \\
> \boldsymbol{u} = \boldsymbol{g} & \text{on } \partial\Omega,
> \end{cases}
> $$
>
> where $\Omega = (0,1)^2$ is the unit square, and $\boldsymbol{g}$ enforces a Poiseuille flow profile with no-slip conditions on the boundary.

The focus of this project is on **preconditioning strategies** for **saddle-point systems**, in particular the **Schur complement method** implemented via PETSc’s `PCFIELDSPLIT`.

## Tools

* **[PETSc](https://petsc.org/)** – Portable, Extensible Toolkit for Scientific Computation
* **[FreeFEM](https://freefem.org/)** – High-level integrated development environment for finite element simulations

## Repository Structure

Each problem has a dedicated folder structured as follows:

```
<PDE_sp_problem>/
├── PETSc/      # C++ implementation using PETSc
└── FreeFEM/    # Scripting implementation using FreeFEM
```

## How to Run

### PETSc

```bash
cd <PDE_sp_problem>/PETSc
cmake -B build
cmake --build build
./build/<PDE_sp_problem> \
  -malloc_dump \
  -ksp_view \
  -ksp_view_final_residual \
  -fieldsplit_0_ksp_converged_reason \
  -fieldsplit_1_ksp_converged_reason
```

Replace `set(PETSC_DIR "~/petsc")` in the CMakeLists.txt with your PETSc installation path.

### FreeFEM

```bash
cd <PDE_sp_problem>/FreeFEM
<PATH_TO_Freefem>/ff-mpirun -np 1 <PDE>.edp \
  -wg -v 0 \
  -malloc_dump \
  -ksp_view \
  -ksp_view_final_residual \
  -fieldsplit_0_ksp_converged_reason \
  -fieldsplit_1_ksp_converged_reason
```

Replace `<PATH_TO_Freefem>` with your FreeFEM installation path, and `<PDE>` with the corresponding `.edp` script. 









