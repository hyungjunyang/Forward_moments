/*
 * File: CYPEqn.h
 * -------------------------------------------------------
 * CYPEqn class, a subclass of Eqn. CYPEqn doesn't own a A, it use 
 * P0Eqn's A by pointing to it. 
 * CYP itself is stored in CYP, size (nNode*nNode)
 *
 * --------------------------------------------------------------------
 *   const double *P0,*dP0dXi[3],*dP0dt -- pointer to P0 and its derivatives
 *   const P0Eqn *P0e                   -- pointer to P0Eqn object
 *
 *   int ic,jc,kc, ijkc                 -- index of reference point
 *   double *CYP                        -- cross-covariance between P1 and Y'
 *
 *  Note: X_old in here, only stores the diagonal part of old CYP
 *
 * --------------------------------------------------------------------
 */

#ifndef _CYPEQN_H
#define _CYPEQN_H
#include <assert.h>
#include <new>
#include "Eqn.h"
#include "P0Eqn.h"

class CYPEqn : public Eqn {
 public:
  CYPEqn( const P0Eqn &P0ee, bool );
  virtual ~CYPEqn();

  // --- setup B, then solve AX=B --------
  virtual void solve(int nWells, const Well *w, const Control &c, 
		     Solver &s, int iter);

  // --- output the CYP at ref Point ---
  virtual void output( const char* fileName, int flag, const Control &c,
		       int slice_index=-1 );
  // [Pipatl]--- output CYP for pressure at ref point and Y field 
  void outputCPY(const char* fileName, int flag, const Control &c,int slice_index=-1 );
  void printCYP();

  // --- Getters ------------------------
  const double* getCYP() const { return CYP; }
  const double* getOldCYPDiag() const { return X_old; }
  const double  getCYP(int i, int ijkc) const 
                { return CYP[i + nNode * ijkc ]; }
	  
 protected:

  // calc B (only internal points)
  virtual void calcRHS( const BType* bType ){;}
  virtual void calcRHS( int & , int & , int &, 
		        const BType *bType, const double *bCond) ;
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
  const double *dP0dXi[3], *dP0dt, *p0; 

  double *CYP;

  void checkBalance();
};

#endif
