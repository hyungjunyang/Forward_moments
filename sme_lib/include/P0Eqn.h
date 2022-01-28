/*
 * File: P0Eqn.h
 * --------------------------------------------------------------------
 * P0Eqn class, a subclass of Eqn. P0Eqn own a A by itself, it need to 
 * allocate and calculate all elements of A. 
 * P0 itself is stored in X, size (nNode)
 *
 * --------------------------------------------------------------------
 *   int nDt              -- number of time steps
 *   double **dP0dXi[3]   -- spacial derivative of P0, record for each t
 *   double **dP0dt       -- time derivative of P0, record for each t
 *   double **P0          -- copy of P0, record for each t
 *
 * --------------------------------------------------------------------
 */
#ifndef _P0EQN_H
#define _P0EQN_H
#include <math.h>

#include "Eqn.h"
#include "Well.h"
#include "Control.h"
#include "Solver.h"


class P0Eqn : public Eqn {
 public:
  P0Eqn( Grid &g, Region *regn, Perm *perm, Fluid *fluid, 
         int nWells, const Well *w, const Control &c, Unit, bool);
  virtual ~P0Eqn();

  // --- setup A and B, then solve AX=B --------
  virtual void solve( int nWells, const Well *w, const Control &c, 
		      Solver &s, int iter );

  void solve(int nWells, const Well *w, const Control &c, 
             Solver &s, int iter, int num_meas, 
	     int* i_pmeas, int* j_pmeas, int* k_pmeas, double*);
  // --- output the P0 ---
  virtual void output( const char* fileName, int flag, const Control &c,
		      int slice_index=-1 );
  
  // --- re-calculate A at the require of other Eqn objects, --- 
  // --- which share the same A with it ---
  void recalcA(const BType *bType, int nWells, const Well *w, double dt);

  // --- Getters ------------------------
  const double* getDP0DXi(int id, int iter) const { return dP0dXi[id][iter]; }
  const double* getDP0Dt(int iter) const { return dP0dt[iter]; }
  const double* getP0() const { return X; }
  const double* getP0(int i) const { return P0[i]; }
 
 protected:

  // initialize  
  virtual void initialize( double X_init );

  // calc A (only internal points)
  virtual void calcA();

  // calc B (only internal points)
  virtual void calcRHS( const BType* bType );

  // calc A, B (boundary points)
  virtual void setBoundaries( bool recalA, const BType *bType,
			      const double *bCond );

  // calc A, B (well blocks) 
  virtual void setWells( bool recalA, int nWells, const Well *w );

  // calc A, B (time-derivatives)
  virtual void calcAccu( bool recalA, double dt );
 
 private:
  int nDt;
  bool debug;

  // for allocating P0's derivatives
  double **dP0dXi[3], **dP0dt, **P0;
  
  // initialize elements in A  
  void initializeA();

  // calc derivatives of P0
  void calcDerivatives( int iter, double dt);

  void calcDerivatives( int iter, double dt, const BType* bType,
			const double* bCond );

  // --- calculate A for boundary points ----
  void aBoundXp( BType bType );
  void aBoundXn( BType bType );
  void aBoundYp( BType bType ); 
  void aBoundYn( BType bType );
  void aBoundZp( BType bType ); 
  void aBoundZn( BType bType );

  // --- calculate B for boundary points ----
  void bBoundXp( BType bType, double bCond );
  void bBoundXn( BType bType, double bCond );
  void bBoundYp( BType bType, double bCond ); 
  void bBoundYn( BType bType, double bCond );
  void bBoundZp( BType bType, double bCond ); 
  void bBoundZn( BType bType, double bCond );
};

#endif
