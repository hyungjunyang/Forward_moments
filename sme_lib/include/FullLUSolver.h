/*
 * File: FullLUSolver.h
 * --------------------------------------------------------------------
 * This is a non-blocked full matrix LU solver
 *  1) Store Full Matrix A(nA, nA), L(nA, nA), U(nA, nA), and A=LU
 *  2) both LU and backSolve are non-blocked, bad cashe usage!!!!!
 *  3) this solver is the slowest and consump most memory, so not prefered
 *
 * --------------------------------------------------------------------
 *   double *A       --  A Matrix, A(nA, nA)
 *   double *L, *U   --  the result of A=LU, L(nA, nA), U(nA, nA)
 *
 * --------------------------------------------------------------------
 * Created    07/24/00          hcao     
 */

#ifndef _FULLLUSOLVER_H
#define _FULLLUSOLVER_H
#include "Solver.h"

class FullLUSolver : public Solver {

 public:
  FullLUSolver(const Eqn &e);
  virtual ~FullLUSolver();

  // --- Multi RHS direct solve -----------------
  virtual void solve(bool doLU, const Eqn &e, int nRHS, double *RHS);
  
 protected:
  virtual void setupA(const Eqn &e);         // build A from Eqn

 private:
  double *A, *L, *U;

  int djlu(double *a, int n, double *l, double *u);             // do LU
  void backSolveL(const double *l, int nRHS, int n, double *x); // solve LY=B
  void backSolveU(const double *u, int nRHS, int n, double *x); // solve UX=Y
};

#endif




