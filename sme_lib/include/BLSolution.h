/*
 * File: BLSolution.h
 * --------------------------------------------------------------------
 * BLSolution moments calculation
 *   1) input travel time mean and variance
 *   2) calculate BLSolution mean and variance
 *
 * --------------------------------------------------------------------
 **Streamline Data***************************************************** 
 *   int nLn, nNd, nPoint -- number of Streamlines, number of point 
 *                           on each streamline, nPoint=nLn*nNd
 *   const double *tauAvg, *tauVar -- travel time mean and variance 
 *   const double *x1C, *etaMnC    -- streamline coord in old Cartesian grid
 *                                    needed on visualization
 * 
 **BLSolution Data*****************************************************
 *   int nS                -- Number of points to discrete BLSolution profile   
 *   double sStar, fpwStar -- front BLSolution and velocity
 *   double viscosR,       -- viscosity ratio (mu_o/mu_w)
 *   double Swc,Sor,       -- end point BLSolution
 *   deltaS = 1-Swc-Sor 
 *   double t              -- require solution at time t
 *   double *satUni        -- deterministic BLSolution solution at time t
 *   double *satAvg        -- BLSolution mean at time t
 *   double *satVar        -- BLSolution variance at time t
 *   double *tau, *s       -- the discreted BLSolution profile(tau vs. sat)  
 *
 * --------------------------------------------------------------------
 * Created    08/03/03         Liyong Li 
 */

#ifndef _BLSolution_H
#define _BLSolution_H
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "FractionFlow.h"
#include <iostream>

class BLSolution {            

 public:
   BLSolution(int, double, double, double);
  ~BLSolution();
  //void get_BL_Table();
  //void get_BL_Table(int, double* tau, double*);
  //void cal_BL_solution(int nPts, double* tau, double* sat);
  void cal_BL_solution();
  double getBLSolution(double);
  double* getTau() {return tau;}
  double getFpwStar() {return fractionFlow->getFpwStar();}
  double getFracFlow(double saturation) {return fractionFlow->fw(saturation);}
  double getFracFlow1(double);
  // --- print() ---
  void print();

 private:
  // --- Fraction Flow Object ---
  double viscosR_, Swc_, Sor_;
  double sStar, fpwStar;
  FractionFlow *fractionFlow;
  
  //BL_1d -- Beckley Leverett Solution for 1D ---
  int nPts_;
  double *tau, *bl_1d, *bl_slope;
  double eps;
  double dtau;
  // --- Private Functions --- 
  void initialize();                                     // initialize
  double root(double a, double eps);                     // root of fpw-a=0
};                                                      
                                                          
#endif

