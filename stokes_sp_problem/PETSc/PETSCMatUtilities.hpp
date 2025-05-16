#ifndef PETSCUTILITIES_H
#define PETSCUTILITIES_H

#include <petsc.h>

/**
 * @brief Assemble the Poisson matrix using sparse AIJ format (CSR: Compressed Sparse Row format).
 * 
 * @param A The PETSC matrix to assemble.
 * @param n The number of grid points.
 */
void assemblePoissonMatrix(Mat &A, int n);

/**
 * @brief Assemble a vector with a specified value.
 * 
 * @param b The PETSC vector to assemble.
 * @param n The number of elements in the vector.
 * @param val The value to set in the vector.
 */
void assembleVector(Vec &b, int n, double val);

/**
 * @brief Assemble the constraint matrix.
 * 
 * @param C The PETSC matrix to assemble.
 * @param n The number of grid points.
 */
void assembleConstraintMatrix(Mat &C, int n);

/**
 * @brief Assemble the saddle point problem using nested matrices and vectors.
 * 
 * @param G The final saddle point matrix (MATNEST type).
 * @param b The combined vector (VecNest type).
 * @param A The Poisson matrix.
 * @param F1 The first vector for the saddle point problem.
 * @param C The constraint matrix.
 * @param D The additional matrix for the saddle point problem.
 * @param F2 The second vector for the saddle point problem.
 * @param n_A The size of matrix A.
 * @param n_D The size of matrix D.
 */
void assembleSaddlePointProblem(Mat &G, Vec &b, Mat &A, Vec &F1, Mat &C, Mat &D, Vec &F2, PetscInt &n_A, PetscInt &n_D);

#endif // PETSCUTILITIES_H
