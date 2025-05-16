#ifndef PETSC_SOLVE_UTILITIES_HPP
#define PETSC_SOLVE_UTILITIES_HPP

#include <petsc.h>

/**
 * @brief Solve a saddle point system using a Schur complement approach.
 * 
 * @param G The saddle point matrix (MATNEST type).
 * @param b The right-hand side vector.
 * @param x The solution vector.
 */
void solveSaddleSystemSchur(Mat G, Vec b, Vec x);

#endif // PETSC_SOLVE_UTILITIES_HPP
