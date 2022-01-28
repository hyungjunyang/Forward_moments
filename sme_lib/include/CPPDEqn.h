/*
 * File: CPPDEqn.h
 * -------------------------------------------------------
 * CPPDEqn class, a subclass of Eqn. CPPDEqn doesn't own a A, it use 
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

#ifndef _CPPDEqn_H
#define _CPPDEqn_H

#include "Eqn.h"
#include "P0Eqn.h"
#include "CYPEqn.h"
#include "CPPEqn.h"
#include "Solver.h"
#include "Control.h"
#include "Well.h"

#include "PressMeasure.h"

class CPPDEqn : public Eqn {
 public:
  //CPPDEqn( const P0Eqn &P0ee, Solver &s, const Control &c);
  CPPDEqn( P0Eqn &P0ee, Solver &s, const Control &c, bool);
  virtual ~CPPDEqn();
  void CPPSolution();
  void variCond();
  void solve();
  void calcCPPMeasbySolveEqn(int nWells, const Well *w, 
                             const Control &c, Solver &s, CYPEqn *CYPee, 
                             CPPEqn *CPPee);
  
  // --- Getters ------------------------
  double* getCPP() { return CPP; }
  const double* getCPP() const { return CPP; }
  const double getCPP(int i, int ijk) const { 
	                     return CPP[i + ijk * nNode]; }
  double getUnCdPStd(int i) {return p_std_UnCd[i];}
  double getCondPStd(int i) {return p_std_Cond[i];}

  void calcUnCdPStd();
  void calcCondPStd();

  void printCYY(double *);
  void printCPP(char*);
  void output();
  void outputObjCondCPP();
  void display(double *);
  
  // --- setup B, then solve AX=B --------
  virtual void solve( int nWells, const Well *w, const Control &c, 
		      Solver &s, int iter ) {;}

  // --- output the CPP at ref Point and diagonal ---
  virtual void output( const char* fileName, int flag, const Control &c,
		       int slice_index=-1 ){;}
 
  void output_PVari( const char* fileName, int flag, int slice_index = -1 ) {;}

 protected:
  
  // calc B (only internal points)
  virtual void calcRHS( int&, const BType* bType, const double *bCond) {;}
  virtual void calcRHS( const BType* bType) {;}
  
  // calc B (well blocks) 
  virtual void setWells( bool recalA, int nWells, const Well *w ) {;}

  // calc B (time-depend parts)
  virtual void calcAccu( bool recalA, double dt ) {;}
 
 private:
  bool debug;
  double small_value;
  int ijk_max;

  // pressure measurements
  int num_p_meas;
  PressMeasure *prsmPtr;
  void addPressMeasureObj(int num);

  //const P0Eqn *P0e;
  P0Eqn *P0e;
  const double *dP0dXi[3], *dP0dt, *p0;

  // the memory required
  double *CYY;
  double *CPP;
  double *p_std_UnCd;
  double *p_std_Cond;
  double *p_std_Meas;

  double *aInv;
  double *b1, *b2, *b3, *b4, *b5, *b6;
  void initialize();
  void readData();
  void calcCYY();
  void YVarInterpolation(double *CYY);
  void calcAInv(Solver &s);
  void calcAmR(int &irow, int &jcol, double *AmB);
  void calcBs(const BType* bType, const double *bCond);
  void calcPermStdBySimplx(int&, double *AmM1, double *AmM3, double *AmO, double *xSolu);
  void simplx(double *amatrix, int *m, int *n, int *mp, int *np, int *m1, 
              int *m2, int *m3, int *icase, int *izrov, int *iposv,
              int *nmax, int *l1, int* mmax, int *l2, int *l3);
};

#endif
