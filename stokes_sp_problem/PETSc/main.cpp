#include "PETSCMatUtilities.hpp"
#include "PETSCSolveUtilities.hpp"

#include "petsc.h"

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    PetscInt n_A;
    PetscInt n_D;

    Mat G, A, C, D;
    Vec b, F1, F2;

    // Assemble the saddle point problem
    assembleSaddlePointProblem(G, b, A, F1, C, D, F2, n_A, n_D);

    // Print the assembly
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, NULL, &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);

    // Solve the linear system
    Vec x;
    VecDuplicate(b, &x);
    solveSaddleSystemSchur(G, b, x);
    
    // Clean up 
    MatDestroy(&G);
    MatDestroy(&A);
    MatDestroy(&C);
    MatDestroy(&D);
    VecDestroy(&b);
    VecDestroy(&F1);
    VecDestroy(&F2);
    VecDestroy(&x);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);
    
    PetscFinalize();
    return 0;
}
