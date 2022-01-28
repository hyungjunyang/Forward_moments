/*
 * File: PVar_Sensitivity.h
 * -------------------------------------------------------
 * PVar_Sensitivity class, a subclass of Eqn. 
 * PVar_Sensitivity doesn't own a A, 
 * it use P0Eqn's A by pointing to it. 
 * PVar_Sensitivity itself is stored in pVar_Sensitivity, size (nNode)
 *
 * --------------------------------------------------------------------
 *   const double *P0, *dP0dXi[3],      -- pointer to P0 and its derivatives
 *   const P0Eqn *P0e                   -- pointer to P0Eqn object
 * --------------------------------------------------------------------
 */

#ifndef _PVar_Sensitivity_H
#define _PVar_Sensitivity_H

#include <new>

#include <assert.h>

#include "Eqn.h"
#include "P0Eqn.h"
#include "CYPEqn.h" 
#include "CPPEqn.h"
#include "Control.h"
#include "Solver.h"

#include "MasterPoint.h"
#include "PressMeasure.h"
#include "KrigWeigh.h"

using namespace std;

class PVar_Sensitivity : public Eqn {
public:
  PVar_Sensitivity(const P0Eqn &P0ee, CYPEqn *CYPe, CPPEqn *CPPee,
                   MasterPoint *mpts, PressMeasure *prsm, KrigWeigh *krig);

  virtual ~PVar_Sensitivity();

  void setStdP();
 
  double getPStd(int i) {return std_PP[i];} 
  double getPVarSensitivity(int msm, int mst) {
	 return pVar_Sensitivity[ msm + mst * num_meas]; }

  // get Functions
  double *getPVarSensitivity() { return pVar_Sensitivity;}

  // --- setup B, then solve AX=B --------
  void setCondPress();
  virtual void solve(int nWells, const Well *w, const Control &c, 
		     Solver &s, int iter);

  void calcCondPStd(int nWells, const Well *w, 
		    const Control &c, Solver &s, int iter );
  void calcCondPStd(int nWells, const Well *w, 
		    const Control &c, Solver &s, double* solu);
  void calcConstrain(int, int, int, double *solu);

  void calcCondPStdByFullEqn(int nWells, const Well *w, 
		             const Control &c, Solver &s, int iter );
  
  void calcSensitivityByFullEqn(int nWells, const Well *w, 
                      const Control &c, Solver &s, int iter);

  void calcSensitity(int nWells, const Well *w, 
                     const Control &c, Solver &s, int iter );
  
  void calcSensitity(int nWells, const Well *w, 
                     const Control &c, Solver &s);
  
  void calcSensiRHS(int&, int & i1_perm, int & j1_perm, int & k1_perm, 
                    const BType *bType, const double *bCond);
  void calcRHS(int & i1_perm, int & j1_perm, int & k1_perm, 
                 const BType *bType, const double *bCond);
  
  virtual void output( const char* fileName, int flag, const Control &c,
		       int slice_index = -1 ) {;} 

  void output();
  void display();
  //void  calcRHO_YP(int nWells, const Well *w, const Control &c, 
  //                  Solver &s, int iter);
  //void solve(const Control &c, Solver &s, int iter) {;}
  //void calcA(const BType *bType, const double *bCond); 
  //void calcRHS( const BType *bType, const double *bCond );

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
  void calcARHO(const BType *bType, const double *bCond, 
		int nWells, const Well *w, int); 
private: 
  bool debug;
  
  const P0Eqn  *P0e;
  const double *dP0dXi[3], *p0; 
  
  CYPEqn *CYPe;
  CPPEqn *CPPe;
	  
  MasterPoint  *mptsPtr;
  PressMeasure *prsmPtr;
  KrigWeigh    *krigPtr;

  int num_meas; 
  double *pVar_Sensitivity;
  double *std_PP, *std_PP_UnCd;
  double *std_YY, *std_YY_UnCd;
  double bigNumber;
  int rho_i, rho_j, rho_k, rho_ijk;
  double *rho_YP;

  void simplx_driver1();
  void simplx_driver2();
  void simplx_search(double upper_bound, double rhs, double *, double *);
		     
  void calcPermStdBySimplx(int& i1_perm, int& j1_perm, int & k1_perm);
  void simplx(double *, int *, int *, int *, int *, 
              int *, int *, int *, int *, int *, int *, 
              int *, int *, int *, int *, int *);
  
  void calcAinv(Solver &s, double *);
  void calcAMs(int i1_perm, int j1_perm, int k1_perm,
               double *am1, double *am2, double *am3, double *am4);
  void calcPermStdBySimplx(double *amR, double *amO, double *); 
};

#endif
