#include "PETSCSolveUtilities.hpp"

void solveSaddleSystemSchur(Mat G, Vec b, Vec x) { // Changed function name
    KSP ksp;
    PC pc;
    PetscInt its;
    PetscReal norm;

    PetscInt m_start, m_end;
    PetscCallVoid(MatGetOwnershipRange(G, &m_start, &m_end));
    PetscPrintf(PETSC_COMM_SELF, " m: %" PetscInt_FMT " n: %" PetscInt_FMT "\n", m_start, m_end);

    // Create linear solver context
    PetscCallVoid(KSPCreate(PETSC_COMM_WORLD, &ksp)); 
    PetscCallVoid(KSPSetOperators(ksp, G, G));

    // OUTER SYSTEM
    // Set up Schur complement preconditioner. Blocks automatically detected since G is MATNEST type.
    PetscCallVoid(KSPGetPC(ksp, &pc));
    PetscCallVoid(PCSetType(pc, PCFIELDSPLIT));
    // Choose Schur decomposition type.
    PetscCallVoid(PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR));
    // Compute full Schur decomposition.
    PetscCallVoid(PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_FULL));
    // Compute Schur complement exactly.
    PetscCallVoid(PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_FULL, NULL));
    // Set the solver for the outer system to PREONLY (the preconditioner is the solver).
    PetscCallVoid(KSPSetType(ksp, KSPPREONLY));

    // INNER SYSTEMS
    // -- FIRST BLOCK
    PetscCallVoid(KSPSetFromOptions(ksp));
    PetscCallVoid(PCSetUp(pc));
    KSP *sub_ksp;
    PetscCallVoid(PCFieldSplitGetSubKSP(pc, NULL, &sub_ksp)); 
    PetscCallVoid(KSPSetType(sub_ksp[0], KSPPREONLY));
    PCSetFromOptions(pc);
    PC sub_pc;
    PetscCallVoid(KSPGetPC(sub_ksp[0], &sub_pc));
    PetscCallVoid(PCSetType(sub_pc, PCLU));

    // -- SECOND BLOCK
    PetscCallVoid(KSPSetType(sub_ksp[1], KSPPREONLY));
    PetscCallVoid(KSPGetPC(sub_ksp[1], &sub_pc));
    // Set the PC for fieldsplit[1] 
    PetscCallVoid(PCSetType(sub_pc, PCLU));

    // SOLVE the linear system
    PetscCallVoid(KSPSolve(ksp, b, x));

    // Get the number of iterations and the norm of the solution
    PetscCallVoid(KSPGetIterationNumber(ksp, &its));
    PetscCallVoid(KSPGetResidualNorm(ksp, &norm));

    // Print the results
    PetscCallVoid(PetscPrintf(PETSC_COMM_WORLD, "Norm of error: %g, Iterations: %" PetscInt_FMT "\n", (double)norm, its)); 

    // Clean up
    PetscCallVoid(KSPDestroy(&ksp));
    PetscCallVoid(PetscFree(sub_ksp));
}
