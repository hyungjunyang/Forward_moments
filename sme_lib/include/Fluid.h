/******************** Property Data ***********************************
 *
 *   double rDen_*  -- reference densities, (w+o) 
 *   double cr      -- rock compressibilities
 *   double cf_w    -- fluid compressibilities (w+o)
 *   double rPoro   -- reference porosity 
 *   double rPres   -- reference pressure
 *   double *poro   -- porosity
 *   double *den_*  -- densities(w+o) for each block 
 *   double *gamma  -- = 9.8 * den 
 *   double *den_cf -- = 9.8*den*cf
 *   double *poro_crf  = poro*(cr+cf)
 *   double *a, *b, *c  -- convinent variables
 *                          (a=9.8*den, b=9.8*den*cf, c=poro*(cr+cf))
 ******************** Saturation Profile Data **************************
 *
 *   double satPType    -- saturation profile type 
 *                         (1:linear; 2:eclipse)
 *   double a0, a1      -- linear : z = a0 + a1*x;
 *                      -- eclipse: (z/a0)^2 + (x/a1)^2 = 1
 *
 * --------------------------------------------------------------------
 */

#ifndef _Fluid_H
#define _Fluid_H

#include <iostream>
#include <fstream>

#include "Common.h"
#include "Grid.h"

class Fluid {
 public:
  Fluid(ifstream &input, Grid*, Unit unit);
 ~Fluid();
  void readData(ifstream &is);
  bool isCompressible() const { return (cr > 0 || cf_o > 0 || cf_w > 0); }
//double getA(int i)          { return a[i]; }
//double getB(int i)          { return b[i]; }
//double getC(int i)          { return c[i]; }
  double getGamma(int i)      { return gamma[i]; }
  double getDen_cf(int i)     { return den_cf[i]; }
  double getPoro_crf(int i)   { return poro_crf[i];}
  double getPoro(int i)       { return poro[i];}
  
  void update(const double *P, bool init);
  
private:
  Unit unit;
  Grid *gridPtr;
  bool debug;
  double cf_o, cf_w, cr;
  double rPoro, rDen_o, rDen_w, rPres;
  double *poro, *den_o, *den_w;
  //double *a,
  double *gamma;
  double *den_cf;
  double *poro_crf;
  //double *b;
  //double *c;
  double a0, a1;
  int satPType;

  void initialize();
  // --- use phase density and saturation profile for the gravity part ---
  void usePhaseInfo_linear();
  void usePhaseInfo_eclipse();
  void convertUnit();
};

#endif
