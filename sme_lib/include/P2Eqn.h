/*
 * File: P2Eqn.h
 * -------------------------------------------------------
 * P2Eqn class, a subclass of Eqn. P2Eqn doesn't own a A, it use 
 * P0Eqn's A by pointing to it. 
 * P2 itself is stored in X, size (nNode)
 *
 * --------------------------------------------------------------------
 *   const double *P0, *dP0dt -- pointer to P0 and its time derivatives
 *   const double *OldCYPDiag -- pointer to old CYP's diagonal part
 *   const double *OldCYP     -- pointer to CYP
 *   const P0Eqn *P0e         -- pointer to P0Eqn object
 *
 *   double *dP2dXi[3]        -- P2's spacial derivatives
 *
 *
 * --------------------------------------------------------------------
 */

#ifndef _P2EQN_H
#define _P2EQN_H
#include "Eqn.h"
#include "P0Eqn.h"
#include "CYPEqn.h"

class P2Eqn : public Eqn {
 public:
  P2Eqn(const P0Eqn &P0ee, const CYPEqn &CYPee, bool);
  virtual ~P2Eqn();

  // --- setup B, then solve AX=B --------
  virtual void solve(int nWells, const Well *w, const Control &c, 
		     Solver &s, int iter);

  // --- output the P0+P2 -----------------
  virtual void output( const char* fileName, int flag, const Control &c,
		       int slice_index=-1 );
  
  // --- Getters ------------------------
  const double* getDP2DXi(int i) const { return dP2dXi[i]; }
  const double* getP2() const { return X; }

 protected:

  // calc B (only internal points
  virtual void calcRHS( const BType* bType ) {;}
  virtual void calcRHS( const BType* bType, const double *bCond );

  // calc B (boundary points)
  virtual void setBoundaries( bool recalA, const BType *bType,
			      const double *bCond ) {;}

  // calc B (well blocks) 
  virtual void setWells( bool recalA, int nWells, const Well *w );

  // calc B (time-depend parts)
  virtual void calcAccu( bool recalA, double dt );

 private:
  bool debug;

  const P0Eqn *P0e;   
  const double *dP0dXi[3];
  const double *dP0dt, *CYP, *OldCYPDiag, *p0;
  double *dP2dXi[3];

  // calc P2's spatial derivatives
  void calcDerivatives(const BType* bType, const double *bCond );

};

#endif

