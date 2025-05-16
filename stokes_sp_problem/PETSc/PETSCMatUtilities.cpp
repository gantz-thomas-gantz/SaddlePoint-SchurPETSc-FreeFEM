#include "PETSCMatUtilities.hpp"

void LoadMatrix(Mat &M, const char *filename) { 
    PetscViewer viewer;
    PetscCallVoid(MatCreate(PETSC_COMM_WORLD, &M));
    char prefix[PETSC_MAX_PATH_LEN];
    PetscCallVoid(PetscSNPrintf(prefix, sizeof(prefix), "%s", filename));
    PetscCallVoid(PetscViewerBinaryOpen(PETSC_COMM_WORLD, prefix, FILE_MODE_READ, &viewer));
    PetscCallVoid(MatLoad(M, viewer));
    PetscCallVoid(PetscViewerDestroy(&viewer));
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

void assembleSaddlePointProblem(Mat &G, Vec &b, Mat &A, Vec &F1, Mat &C, Mat &D, Vec &F2, PetscInt &n_A, PetscInt &n_D) {
    Mat Ct;
    PetscViewer viewer;

    LoadMatrix(A, "../data/B00.dat");
    LoadMatrix(C, "../data/B10.dat");
    LoadMatrix(D, "../data/B11.dat");

    // Create transpose of C
    PetscCallVoid(MatTranspose(C, MAT_INITIAL_MATRIX, &Ct));

    // Create nested matrix
    Mat nestedMats[4] = {A, Ct, C, D};
    PetscCallVoid(MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, &nestedMats[0], &G));
    PetscCallVoid(MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY));
    PetscCallVoid(MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY));

    // Combine vectors into a nested vector
    PetscCallVoid(MatGetSize(A, &n_A, NULL));
    PetscCallVoid(MatGetSize(D, &n_D, NULL));
    assembleVector(F1, n_A, 0);
    assembleVector(F2, n_D, 0);
    Vec nestedVecs[2] = {F1, F2};
    PetscCallVoid(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, nestedVecs, &b));
    
    PetscCallVoid(MatDestroy(&Ct));
}
