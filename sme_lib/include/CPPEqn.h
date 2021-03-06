/*
 * File: CPPEqn.h
 * -------------------------------------------------------
 * CPPEqn class, a subclass of Eqn. CPPEqn doesn't own a A, it use 
 * P0Eqn's A by pointing to it. 
 * CPP itself is stored in CPP, size (nNode*nNode)
 *
 * --------------------------------------------------------------------
 *   const double *P0,*dP0dXi[],*dP0dt --  pointer to P0 and its derivatives
 *   const double *CYP                 -- pointer to CYP
 *   const P0Eqn  *P0e                 -- pointer to P0Eqn object
 *
 *   int ic,jc,kc, ijkc                -- index of reference point
 *   double *CPP                       -- covariance of P1
 *
 * --------------------------------------------------------------------
 */

#ifndef _CPPEQN_H
#define _CPPEQN_H
#include "Eqn.h"
#include "P0Eqn.h"
#include "CYPEqn.h"

class CPPEqn : public Eqn {
 public:
  CPPEqn( const P0Eqn &P0ee, const CYPEqn &CYPee, bool );
  virtual ~CPPEqn();

  // --- setup B, then solve AX=B --------
  virtual void solve( int nWells, const Well *w, const Control &c, 
		      Solver &s, int iter );
  
  // --- output the CPP at ref Point and diagonal ---
  virtual void output( const char* fileName, int flag, const Control &c,
		       int slice_index=-1 );
 
  void output_PVari( const char* fileName, int flag, int slice_index = -1 );
 
  // --- Getters ------------------------
  double* getCPP() { return CPP; }
  const double* getCPP() const { return CPP; }
  const double getCPP(int i, int ijk) const {
	                     return CPP[i + ijk * nNode]; }
  double getUnCdPStd(int i) {return p_std_UnCd[i];}
  double getCondPStd(int i) {return p_std_Cond[i];}
  double getUnCdRho(int i)  {return rho_YP_UnCd[i];}
  double getCondRho(int i)  {return rho_YP_Cond[i];}
  void calcUnCdPStd();
  void calcCondPStd();
  void printCPP();
  void printRho();
 
 protected:
  
  // calc B (only internal points)
  virtual void calcRHS( int&, const BType* bType, const double *bCond);
  virtual void calcRHS( const BType* bType) {;}
  
  // calc B (well blocks) 
  virtual void setWells( bool recalA, int nWells, const Well *w ) {;}

  // calc B (time-depend parts)
  virtual void calcAccu( bool recalA, double dt );
 
 private:
  bool debug;
  const P0Eqn *P0e;
  const double *dP0dXi[3], *dP0dt, *CYP, *p0;

  double *CPP;
  double *p_std_UnCd;
  double *p_std_Cond;
  double *rho_YP_UnCd;
  double *rho_YP_Cond;

  void checkBalance();
};

#endif
