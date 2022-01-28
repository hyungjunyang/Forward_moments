/*
 * File: LapackBandSolver.h
 * -------------------------------------------------------
 * This is a blocked band matrix LU solver from LAPACK
 *  1) Store band Matrix A in special format, storage A(3*bandsize+1, nA) 
 *  2) After LU, L and U are stored inside A, 
 *       please refer to Lapack menu (dgbtrf) for the detail storage format!!!
 *  3) both LU and backSolve are blocked, better cashe usage.
 *  4) this solver is much faster than Full Matrix Solver due to the band
 *     format of the matrix, it also consump much less memory. 
 *  5) The performance decreses with bandsize increasing, after bandsize 
 *     increases to about a quarter of A's size(nA), its performance is worse 
 *     than full matrix Lapack solver. Stop Using it!!!!!!!!!!!!
 *  6) Lapack subroutine dgbtrf(LU) and dgbtrs(back solve) are used.
 *  7) This code implement the bandTransform process to minimize bandsize.
 *           1D: bandsize=1
 *           xy: badnsize=nx, so best when nx<=ny
 *           xz: bandsize=min(nx, nz), always best  (bandTransform if nx>nz)
 *           3D: bandsize=ny*nz if nx>nz (best when nx>=ny) (bandTransform)
 *               bandsize=ny*nx if nx<=nz  (no bandTransform)
 *      In summary: bandTransform will switch I+- bands with K+- bands 
 *      if required(nz>1&&nx>nz).
 *
 * --------------------------------------------------------------------
 *   int bandSize    --  A's band size (band width)
 *   int aHeight     --  A's number of rows = 3*bandSize+1
 *   bool bandSwitch --  decide do bandTransform or not
 *   double *A       --  A Matrix, A(aHeight, nA)
 *   int *ipiv       --  pivoting information produced by LU, 
 *                                     and used by backSolve
 *
 * --------------------------------------------------------------------
 * Created    07/25/00          hcao     
 */

#ifndef _LAPACKBANDSOLVER_H
#define _LAPACKBANDSOLVER_H

#include "Solver.h"

class LapackBandSolver : public Solver {

 public:
  LapackBandSolver(const Eqn &e);
  virtual ~LapackBandSolver();

  // --- Multi RHS direct solve -----------------
  virtual void solve(bool doLU, const Eqn &e, int nRHS, double *RHS);
   
 protected:
  virtual void setupA(const Eqn &e);     // build A from Eqn (no band transform) 
  void setupA_Transform(const Eqn &e);   // build A from Eqn (band transform)

 private:
  int  bandSize, aHeight;
  bool bandSwitch;
  int *ipiv;
  double *A;

  // --- Wrappers to Lapack Fortran subroutines (dgbtrf, dgbtrs) ------
  int DGBTRF_C(int m, int n, int kl, int ku, double *a, int lda, int *ipiv);
  int DGBTRS_C(char tran, int n, int kl, int ku, int nRHS, 
	       double *a, int lda, int *ipiv, double *B, int ldb);
};

#endif




