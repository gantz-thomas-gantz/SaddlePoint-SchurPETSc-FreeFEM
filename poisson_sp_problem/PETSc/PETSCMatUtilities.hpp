#ifndef PETSCUTILITIES_H
#define PETSCUTILITIES_H

// Include PETSC header
#include <petsc.h>

/**
 * @brief Assemble the Poisson matrix using sparse AIJ format (CSR: Compressed Sparse Row format).
 * 
 * @param A The PETSC matrix to assemble.
 * @param n The number of grid points.
 */
void assemblePoissonMatrix(Mat &A, int n);

/**
 * @brief Assemble a vector.
 * 
 * @param b The PETSC vector to assemble.
 * @param n The number of elements in the vector.
 */
void assembleVector(Vec &b, int n, int val);

/**
 * @brief Assemble the constraint matrix.
 * 
 * @param C The PETSC matrix to assemble.
 * @param n The number of grid points.
 */
void assembleConstraintMatrix(Mat &C, int n);

/**
 * @brief Assemble the saddle point problem.
 * 
 * @param G The final saddle point matrix.
 * @param b The combined vector.
 * @param A The Poisson matrix.
 * @param F The Poisson vector.
 * @param C The constraint matrix.
 * @param c The constraint vector.
 * @param n The number of grid points.
 * @param c_val The constraint value.
 */
void assembleSaddlePointProblem(Mat &G, Vec &b, Mat &A, Vec &F, Mat &C, Vec &c, int n, double c_val);

#endif // PETSCUTILITIES_H
