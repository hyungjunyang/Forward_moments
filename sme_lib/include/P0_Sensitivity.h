/*
 * File: P0_Sensitivity.h
 * -------------------------------------------------------
 * P0_Sensitivity class, a subclass of Eqn. 
 * P0_Sensitivity doesn't own a A, it use 
 * P0Eqn's A by pointing to it. 
 * P0_Sensitivity itself is stored in CYP, size (nNode*nNode)
 *
 * --------------------------------------------------------------------
 *   const double *P0,*dP0dXi[3],*dP0dt -- pointer to P0 and its derivatives
 *   const P0Eqn *P0e                   -- pointer to P0Eqn object
 *
 *   int ic,jc,kc, ijkc                 -- index of reference point
 *
 *  Note: X_old in here, only stores the diagonal part of old CYP
 *
 * --------------------------------------------------------------------
 */

#ifndef _P0_Sensitivity_H
#define _P0_Sensitivity_H
#include <assert.h>
#include <new>

#include "Eqn.h"
#include "P0Eqn.h"
#include "Control.h"
#include "Solver.h"

#include "MasterPoint.h"
#include "PressMeasure.h"
#include "KrigWeigh.h"

class P0_Sensitivity : public Eqn {
public:
  P0_Sensitivity(const P0Eqn &P0ee, MasterPoint *mpts,
                 PressMeasure *prsm, KrigWeigh *krig);

  virtual ~P0_Sensitivity();

  double getP0Sensitivity(int msm, int mst) {
	  return p0_Sensitivity[ msm + mst * num_meas]; }

  // get Functions
  double *getP0Sensitivity() { return p0_Sensitivity;}

  // --- setup B, then solve AX=B --------
  virtual void solve(int nWells, const Well *w, const Control &c, 
		     Solver &s, int iter) {;}
  void solve(const Control &c, Solver &s);
    
  virtual void calcRHS(int& master, const BType *bType, const double *bCond);

  virtual void output( const char* fileName, int flag, const Control &c,
		       int slice_index = -1 ) {;} 

  void output();
  void display();

 protected:

  // calc B (only internal points)
  virtual void calcRHS( const BType* bType ){;}
  
  // calc B (boundary points)
  virtual void setBoundaries( bool recalA, const BType *bType,
			      const double *bCond ) {;}

  // calc A, B (well blocks) 
  virtual void setWells( bool recalA, int nWells, const Well *w ){;}

  // calc A, B (time-depend parts)
  virtual void calcAccu( bool recalA, double dt );

private: 
  bool debug;

  const P0Eqn *P0e; 
  const double *dP0dXi[3], *p0; 

  MasterPoint *mptsPtr; 
  PressMeasure *prsmPtr;
  KrigWeigh   *krigPtr;

  int num_meas;
  double *p0_Sensitivity;
};

#endif
