#include <iostream>
#include "PETSCMatUtilities.hpp"
#include "PETSCSolveUtilities.hpp"

#include "petsc.h"

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    int n = 10;         // Number of grid points
    double h = 1./n;    // Grid spacing
    double c_val = 1.0; // Average over volume constraint value

    Mat G, A, C;
    Vec b, F, c;

    // Assemble the saddle point problem
    assembleSaddlePointProblem(G, b, A, F, C, c, n, c_val);

    // Print the assemblage
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, NULL, &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);

    // std::cout << "Poisson Matrix A:" << std::endl;
    // MatView(A, viewer);

    // std::cout << "Saddle Point Matrix G:" << std::endl;
    // MatView(G, viewer);

    // std::cout << "Combined Vector b:" << std::endl;
    // VecView(b, viewer);

    // Solve the linear system
    Vec x;
    VecDuplicate(b, &x);
    solveSaddleSystemSchur(G, b, x);
    
    // std::cout << "Solution x:" << std::endl;
    // VecView(x, viewer);
    
    // Clean up 
    MatDestroy(&G);
    MatDestroy(&A);
    MatDestroy(&C);
    VecDestroy(&b);
    VecDestroy(&F);
    VecDestroy(&c);
    VecDestroy(&x);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);
    
    PetscFinalize();
    return 0;
}
