/*
 * File: Solver.h
 * --------------------------------------------------------------------
 * Base class for all type of the direct solvers used in this simultor
 *   Use direct method to solve AX=B, could be multi-RHS B.
 *   1) Eqn class provides the "raw" data for A, B and X
 *   1) A is rebuild (copy one by one) according to specific format 
 *      required by each solver, the diag and offdiag bands in Eqn class 
 *      are passed in to setup A. 
 *   2) X and B share the same space. 
 *      Before solve, it store B, after solver, it store X.
 *      X|B directly use data from Eqn class (pass pointer), no rebuild
 *   3) A is stored in subclasses (individual solvers), 
 *      which know the format of A
 *
 * --------------------------------------------------------------------
 *   int nx, ny, nz, nNode -- grid dimension, nNode=nx*ny*nz
 *   int nW, nA            -- number of CONST_RATE multi-block penetrated 
 *                            wells, which will add a new line to Jacobi.
 *                            nA = nNode+nW, the size if A, 
 *                            A(nA, nA), if A is a full matrix
 *
 * --------------------------------------------------------------------
 * Created    06/22/00          hcao     
 */

#ifndef _SOLVER_H
#define _SOLVER_H
#include "Common.h"

class Eqn;
class Solver {

public:
  Solver(const Eqn& e);   
  virtual ~Solver();

  // --- Multi RHS direct solve -----------------
  virtual void solve(bool doLU, const Eqn &e, int nRHS, double *RHS) = 0;
  
protected:
  bool debug;
  int nx, ny, nz, nNode, nA, nW;
  int nxPerm, nyPerm, nzPerm, nNodePerm;
  
  virtual void setupA(const Eqn &e) = 0;     // build A from Eqn
};

#endif
