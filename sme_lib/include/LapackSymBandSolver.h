/*
 * File: LapackSymBandSolver.h
 * -------------------------------------------------------
 * This is a blocked symmetric band matrix LU solver from LAPACK
 *  1) Store band Matrix A in special format, storage A(bandsize+1, nA) 
 *  2) After LU, L and U are stored inside A, 
 *       please refer to Lapack menu (dpbtrf) for the detail storage format!!!
 *  3) both LU and backSolve are blocked, a little better cashe usage.
 *  4) this solver is much faster than Full Matrix Solver due to the band
 *     format of the matrix, and slightly faster than band Matrix Solver 
 *     due to the symmetric condition, it consumpts least memory. 
 *  5) The performance decreses with bandsize increasing, after bandsize 
 *     increases to about a quarter of A's size(nA), its performance is worse 
 *     than full matrix Lapack solver. Stop Using it!!!!!!!!!!!!
 *  6) Lapack subroutine dpbtrf(Ches..) and dpbtrs(back solve) are used.
 *  7) Ths solver only works for positive defined symmetric band matrix.
 *  8) This code also implements the bandTransform process to minimize 
 *     bandsize, detail refer to LapackBandSolver.h
 *
 * --------------------------------------------------------------------
 *   int bandSize    --  A's band size (band width)
 *   int aHeight     --  A's number of rows = bandSize+1
 *   bool bandSwitch --  decide do bandTransform or not
 *   double *A       --  A Matrix, A(aHeight, nA)
 *        no need for pivoting due to positive defined condition...
 *
 * --------------------------------------------------------------------
 * Created    07/27/00          hcao     
 */

#ifndef _LAPACKSYMBANDSOLVER_H
#define _LAPACKSYMBANDSOLVER_H

#include "Solver.h"

class LapackSymBandSolver : public Solver {

 public:
  LapackSymBandSolver(const Eqn &e);
  virtual ~LapackSymBandSolver();

  // --- Multi RHS direct solve -----------------
  virtual void solve(bool doLU, const Eqn &e, int nRHS, double *RHS);
 
 protected:
  virtual void setupA(const Eqn &e);  // build A from Eqn (no band transform) 
  void setupA_Transform(const Eqn &e);// build A from Eqn (band transform)

 private:
  int bandSize, aHeight;
  bool bandSwitch;
  double *A;

  // --- Wrappers to Lapack Fortran subroutines (dpbtrf, dpbtrs) ------
  int DPBTRF_C(char uplo, int n, int kd, double *a, int lda);
  int DPBTRS_C(char uplo, int n, int kd, int nRHS, 
	       double *a, int lda, double *B, int ldb);
};

#endif




