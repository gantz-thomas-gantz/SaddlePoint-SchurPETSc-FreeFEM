#include "PETSCMatUtilities.hpp"

void assemblePoissonMatrix(Mat &A, int n) {
    double h = 1.0/n;
    PetscCallVoid(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n, 3, NULL, 3, NULL, &A));

    PetscInt m_start, m_end;
    PetscCallVoid(MatGetOwnershipRange(A, &m_start, &m_end));
    for (int i = m_start; i < m_end; ++i) {
        if (i > 0) {
            PetscCallVoid(MatSetValue(A, i, i-1, -1.0/h, INSERT_VALUES));
        }
        PetscCallVoid(MatSetValue(A, i, i, 2.0/h, INSERT_VALUES));
        if (i < n-1) {
            PetscCallVoid(MatSetValue(A, i, i+1, -1.0/h, INSERT_VALUES));
        }
    }
    PetscCallVoid(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCallVoid(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
}

void assembleVector(Vec &b, int n, double val) {
    PetscCallVoid(VecCreate(PETSC_COMM_WORLD, &b));
    PetscCallVoid(VecSetSizes(b, PETSC_DECIDE, n));
    PetscCallVoid(VecSetFromOptions(b));
    PetscInt m_start, m_end;
    PetscCallVoid(VecGetOwnershipRange(b, &m_start, &m_end));
    for (int i = m_start; i < m_end; ++i) {
        PetscCallVoid(VecSetValue(b, i, val, INSERT_VALUES));
    }
    PetscCallVoid(VecAssemblyBegin(b));
    PetscCallVoid(VecAssemblyEnd(b));
}

void assembleConstraintMatrix(Mat &C, int n) {
    double h = 1.0/n;
    PetscCallVoid(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, 1, 1, NULL, 1, NULL, &C));
    PetscInt m_start, m_end;
    PetscCallVoid(MatGetOwnershipRange(C, &m_start, &m_end));
    for (int i = m_start; i < m_end; ++i) {
        PetscCallVoid(MatSetValue(C, i, 0, h, INSERT_VALUES));
    }
    PetscCallVoid(MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY));
    PetscCallVoid(MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY));
}

void assembleSaddlePointProblem(Mat &G, Vec &b, Mat &A, Vec &F, Mat &C, Vec &c, int n, double c_val) {
    Mat Ct, ZERO;

    double h = 1.0/n;
    // Assemble submatrices and subvectors
    assemblePoissonMatrix(A, n);
    assembleVector(F, n, h);
    assembleVector(c, 1, c_val);

    // Assemble constraint matrix
    assembleConstraintMatrix(C, n);

    // Create transpose of C
    PetscCallVoid(MatTranspose(C, MAT_INITIAL_MATRIX, &Ct));

    // Create nested matrix
    Mat nestedMats[4] = {A, C, Ct, NULL};
    PetscCallVoid(MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, &nestedMats[0], &G));
    PetscCallVoid(MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY));
    PetscCallVoid(MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY));

    // Combine vectors into a nested vector
    Vec nestedVecs[2] = {F, c};
    PetscCallVoid(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, nestedVecs, &b));
    
    PetscCallVoid(MatDestroy(&Ct));
}
