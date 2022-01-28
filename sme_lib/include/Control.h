/*
 * File: Control.h
 * --------------------------------------------------------------------
 * Control class defines the parameters and options used in this simultor
 * Including: boundary condition, initial condition, 
 *            saturation profile (front), solver choice, time steps,
 *
 * --------------------------------------------------------------------
 **Boundary Condition: (3D has 6 surface, ordered in X+, X-, Y+, Y-, Z+, Z-) 
 *   BType bType[6]  -- boundary type (1:CONST_RATE, 2:CONST_PRES)  
 *   double bCond[6] -- boundary cond (pres for const pres boundary,
 *                                     rate for const rate boundary)
 *      Note: rate is for each cell on the boundary, not the total
 *
 **Initial Condition: (const pressure everywhere)**********************
 *   int pInit      -- One pressure assign to every PO,
 *                     the initial value for CYP, CPP, P2 are all zero
 *
 **Saturation Profile: (only used for gravity in the z direction)******
 *   double a0, a1  -- z=a0+a1*x (line) for linear flow,
 *                     (z/a0)^2+(x/a1)^2=1 (eclipse) for 5-spot   
 *
 **Solver Choice: (currently 4 choice)*********************************
 *   SolverType solverChoice  -- 
 *   (1:FULLLU, 2:LAPACK_BAND, 3:LAPACK_FULL, 4:LAPACK_SYM_BAND, 5 CLUB_SOLVER)
 *                    
 **Time Steps:********************************************************* 
 *   int nDt        -- number of time steps
 *   double *dt     -- dt for each time step
 *     Note: For incompressible flow, nDt will be set to "1" automatically
 *                                    and dt[0] will be "0.0"
 *
 **Reference point for output:*****************************************
 *   int icc, jcc, kcc  -- index of the reference point
 *
 * --------------------------------------------------------------------
**/

#ifndef _CONTROL_H
#define _CONTROL_H
#include <fstream>
#include <stdlib.h>
#include "Common.h"
#include "Grid.h"
#include "Point.h"
#include <iostream>

class Control {

public:
  // --- Constructor and destructor
  Control(ifstream &is, Grid*, Unit unit);
  ~Control();

  // --- Get number of time steps
  int     getI_Ref()   const { return icc;}
  int     getJ_Ref()   const { return jcc;}
  int     getK_Ref()   const { return kcc;}
  int     getNDt()     const { return nDt;}
  double  getP_Init()  const { return pInit;}
  double  getDt(int i) const { return dt[i];}
  double* getBCond()   const { return bCond;}
  BType*  getBType()   const { return bType;}
  
  double getSwc()     {return Swc;}
  double getSor()     {return Sor;}
  double getViscosR() {return viscosR;}
  double getProductionTime()     {return time_production;}
  int    getProductionTimeStep() {return num_time_step;}
  int    getFlagPressCond()      {return flag_press_condition; } 
  int    getFlagSatCond()		 {return flag_sat_condition;} //Pipatl
  int	 getFlagProdCond()		 {return flag_prod_condition;} //Pipatl
  bool   isInjcBound( int i )  { return bInjc[i] == 1 ? true : false; } //hjyang
  bool   isProdBound( int i )  { return bProd[i] == 1 ? true : false; } //hjyang

  double  getBoundLength( int i ) { return bLength[i]; } 
  Point* getBoundPts( int i ); 
  // --- Get linear solver type
  SolverType getSolverChoice() const { return solverChoice; }

  // --- Set steady state problem for incompressible flow
  void setSteady() { nDt = 1; dt[0] = 0.0; }
 
  // --- output the control parameters
  friend ostream &operator<< (ostream &os, Control &con);
 
private:
  Grid *gridPtr; 
  bool debug;
  int nDt, icc, jcc, kcc;
  double pInit, *dt;
 
  // Production
  double Swc, Sor;
  double viscos_w, viscos_o, viscosR;
  double time_production;
  int    num_time_step;

  // Pressure Conditioning
  int flag_press_condition;

  // Saturation Conditioning (Pipatl)
  int flag_sat_condition;

  // Production Conditioning (Pipatl)
  int flag_prod_condition;

  //double bCond[6],
  //BType  bType[6];
  double *bCond;
  double *bLength; 
  Point *bPtSetXn; 
  Point *bPtSetXp; 
  Point *bPtSetYn;
  Point *bPtSetYp;
  
  double *bInjc; 
  double *bProd;

  BType  *bType; 
  SolverType solverChoice;

  // --- read the input file
  void readData(ifstream &is);
  // --- convert the field unit to metric unit
  void convertUnit();
  // --- calculate the SL points at boundary 
  void calcBoundPoint(); 
};

#endif
