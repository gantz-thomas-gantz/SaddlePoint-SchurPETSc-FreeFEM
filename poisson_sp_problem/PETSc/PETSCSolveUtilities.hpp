#ifndef PETSC_SOLVE_UTILITIES_HPP
#define PETSC_SOLVE_UTILITIES_HPP

#include <petsc.h>

void solveSaddleSystemSchur(Mat G, Vec b, Vec x);

#endif // PETSC_SOLVE_UTILITIES_HPP
