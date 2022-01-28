/*
 * File: Velocity.h
 * --------------------------------------------------------------------
 * This is velocity class, calculate velocity moments(v0, v2, Cvv)  
 *
 * --------------------------------------------------------------------
 *   Unit unit             -- flag for unit (1:FIELD; 2:METRIC)
 *   int nx, ny, nz, nNode -- grid dimension, nNode=nx*ny*nz
 *   const double *dP0dXi[3],*dP2dXi[3] -- pointer to P0 and P2's derivatives
 *   const double *CYP                  -- pointer to CYP
 *   const double *CYP                  -- pointer to CYP
 *   const P0Eqn  *P0e                  -- pointer to P0Eqn object
 *   const Grid   &g                    -- refernce to Grid object for easy access
 *
 *   double *v0[3]         -- velocity mean
 *   double *v2[3]         -- second order velocity 
 *   double *Cv            -- one (cross)-covariance between v1
 *
 * --------------------------------------------------------------------
 */

#ifndef _VELOCITY_H
#define _VELOCITY_H
#include <stdio.h>
#include <stdlib.h>

#include "Grid.h"
#include "Region.h"
#include "Perm.h"
#include "Fluid.h"
#include "Control.h"
#include "P0Eqn.h"
#include "CYPEqn.h"
#include "CPPEqn.h"
#include "P2Eqn.h"
  
class Grid;
class Control;

class Velocity {
 public:
  Velocity(Grid*, Region *, Perm *, Fluid *, 
           const P0Eqn &P0e, const CYPEqn &CYPe, const CPPEqn &CPPe, 
           const P2Eqn &P2e, Unit unit, char *Dir, bool);
  ~Velocity();
  
  // --- solve for velocity mements ---------------
  void solve(const Control &con, int iter); 
  
  double* getV0(int id) { return v0[id]; }
  double* getCv1v1() {return Cv1v1;}
  double* getCv2v2() {return Cv2v2;}
  double* getCv3v3() {return Cv2v2;}
  double* getCv1v2() {return Cv1v2;}
  double* getCv1v3() {return Cv1v1;}
  double* getCv2v3() {return Cv2v2;}

  void check(); 
  void check3d();
  void convertUnit();
  void writeCv(const Control &c);
  
 private:
  bool debug;

  char *directory;  
  Unit unit;
  
  Grid    *gridPtr;
  Region  *regnPtr;
  Perm    *permPtr;
  Fluid   *fluidPtr;

  const P0Eqn *P0e;
  const double *dP0dXi[3], *CYP, *CPP;
  
  int nx, ny, nz, nNode;
  double *v0[3]; 
  double *Cv1v1, *Cv2v2, *Cv3v3;
  double *Cv1v2, *Cv1v3, *Cv2v3;
  
  void initialize();
  void calcVs();
  void calcCViVj( double *, double *, double *);

};

#endif
