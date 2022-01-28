/*
 * File: LapackFullSolver.h
 * -------------------------------------------------------------------
 *  This is a blocked full matrix LU solver from LAPACK
 *  1) Store Full Matrix A(nA, nA), 
 *  2) After LU, L is at the lower tri of A, U is at the upper tri of A
 *  3) both LU and backSolve are blocked, best cashe usage.
 *  4) this solver is faster than FullLUSolver due to cashe usage,
 *     but it still consumpts a lot memory, 
 *  5) when A is not in band format (due to well), it is the best choice.
 *  6) Lapack subroutine dgetrf(LU) and dgetrs(back solve) are used.
 *
 * --------------------------------------------------------------------
 *   double *A       --  A Matrix, A(nA, nA)
 *   int *ipiv       --  pivoting information produced by LU, 
 *                                     and used by backSolve
 *
 * --------------------------------------------------------------------
 * Created    07/27/00          hcao     
 */

#ifndef _LAPACKFULLSOLVER_H
#define _LAPACKFULLSOLVER_H

#include "Solver.h"

class LapackFullSolver : public Solver {

 public:
  LapackFullSolver(const Eqn &e);
  virtual ~LapackFullSolver();

  // --- Multi RHS direct solve -----------------
  virtual void solve(bool doLU, const Eqn &e, int nRHS, double *RHS);

 protected:
  virtual void setupA(const Eqn &e);        // build A from Eqn
  
 private:
  int *ipiv;
  double *A;

  // --- Wrappers to Lapack Fortran subroutines (dgetrf, dgetrs) ------
  int DGETRF_C(int m, int n, double *a, int lda, int *ipiv);
  int DGETRS_C(char tran, int n, int nRHS, double *a, int lda,
	       int *ipiv, double *b, int ldb);
};

#endif




