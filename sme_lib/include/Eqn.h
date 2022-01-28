/*
 * File: Eqn.h
 * --------------------------------------------------------------------
 * This is the base Eqn class, the parent for all individual AX=B type 
 * Eqn classes, currently including P0Eqn, CYPEqn, CPPEqn and P2Eqn
 *
 * Note: 
 *   All terms in Matrix A are located in subclasses according to the need.
 *   If several subclasses share the same A, then only one A is located, 
 *   and each of these classes keeps a pointer to A's elements 
 *       (including diagA, offDA*, WW, WR, RW, wEqnInd) 
 *
 * --------------------------------------------------------------------
 *   const Grid &g         -- refernce to Grid object for easy access
 *   Unit unit             -- flag for unit (1:FIELD; 2:METRIC)
 *   int nx, ny, nz, nNode -- grid dimension, nNode=nx*ny*nz
 *
 *   int nW, nA            -- number of CONST_RATE multi-block penetrated 
 *                            wells, which will add a new line to Jacobi.
 *                            nA = nNode+nW, the size if A, 
 *   double *RHS, *X, *X_old  -- storage for B, X and X at old time step
 *   double *diag, *offD*  -- diagonal and off-diagonals for reservoir 
 *                            part of Matrix A
 *   int *wEqnInd          -- mapping from well's order in A to well's 
 *                            order in original wells 
 *   double *WW, *WR, *RW  -- well parts of A 
 *                            (well to well, well to res, res to well)
 *   const Grid &g         -- refernce to Grid object for easy access
 *
 * --------------------------------------------------------------------
 * Created    06/22/00          hcao     
 */
#ifndef _EQN_H
#define _EQN_H

#include <iomanip>
#include <math.h>

#include "Grid.h"
#include "Region.h"
#include "Perm.h"
#include "Fluid.h"
#include "Common.h"
#include "Control.h"
#include "Well.h"
#include "Solver.h"

class Control;
class Well;
class Solver;

class Eqn {
public:
  Eqn( Grid &grid, Region *regn, Perm *perm, Fluid *fluid );        
  // for subclasses which own an A
  Eqn(const Eqn &e);            // for subclasses which doesn't own an A
  virtual ~Eqn();

  // --- setup A and B, then solve AX=B --------
  virtual void solve(int nWells, const Well *w, const Control &c, 
		     Solver &s, int iter) = 0;
  // --- output the results ---
  virtual void output( const char* fileName, int flag, const Control &c,
		       int slice_index = -1) = 0;
  
  friend ostream &operator << (ostream &os, Eqn &eqn);

  friend class Solver;
  friend class FullLUSolver;
  friend class LapackBandSolver;
  friend class LapackFullSolver;
  friend class LapackSymBandSolver;
  friend class ClubSolver;
  
 protected:
  Unit unit;
  Grid   &g;
  Region *regnPtr;
  Perm   *permPtr;
  Fluid  *fluidPtr;
  int nx, ny, nz, nNode, nW, nA, *wEqnInd; 
  int nxPerm, nyPerm, nzPerm, nNodePerm;
  double *RHS, *X, *X_old;                   // (nA, 1)

  // reservoir parts of A 
  double *diagA;
  double *offDAXp, *offDAYp, *offDAZp;
  double *offDAXn, *offDAYn, *offDAZn; 

  // (nW,nW), (nNode,nW), (nW,nNode)                                          
  // well parts of A
  double *WW, *WR, *RW;  

  // initialize and initial value 
  virtual void initialize( double X_init );

  // calc A (only internal points)
  virtual void calcA() { };

  // calc B (only internal points)
  virtual void calcRHS( const BType* bType) = 0;

  // calc A, B (boundary points)
  virtual void setBoundaries( bool recalA, const BType *bType,
			      const double *bCond );

  // calc A, B (well blocks)
  virtual void setWells( bool recalA, int nWells, const Well *w) = 0;

  // calc A, B (time-depend parts)
  virtual void calcAccu( bool recalA, double dt) = 0;
};

#endif
