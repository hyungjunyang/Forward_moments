/*
 * File: Well.h
 * --------------------------------------------------------------------
 * Well class only deal with the well data input and output, and the basic 
 * calculations, such as q=Tw(Pw-P), Pw=P+q/Tw... 
 * All of the Jacobi and RHS calculations are in Eqn class.
 * Simple well treatment, and no WellManager is defined
 *
 * --------------------------------------------------------------------
 *  WConType wConType -- well control type (1:RATE; 2:PRES) 
 *  int wLens         --  number of penetrations
 *  int *wLocIJK[3]   -- i,j,k index for each well block
 *  double *wTrans    -- well transmissiblity for each well block
 *  double wStr       -- pres for const pres well, rate for const rate well
 *  
 * --------------------------------------------------------------------
 */
#ifndef _WELL_H
#define _WELL_H
#include <fstream>
#include <stdlib.h>
#include "Common.h"
#include <iostream>

class Well {
 public:
   Well();
  ~Well();
   
  void readData(ifstream &is, Unit unit); 
  void update();         // update Pw for RATE, and q for PRES
 
  WConType getWellType()            const { return wConType; } 
  int      getWellLen()             const { return wLens; }
  double   getWellForce()           const { return wStr;}
  double   getWellTrans(int i)      const { return wTrans[i];} 
  int      getWellIJK(int i, int j) const { return wLocIJK[i][j];}
  double   getWellXYZ(int i, int j) { return wLocXYZ[i][j];}
  int      getNumPoints()   {return npts;}
  double   getDegree(int i) {return angle[i];}
  double   getRadius()      {return radius; }
  bool isProd() { return isProdWell; }
  bool isInjc() { return isInjcWell; }

  friend ostream &operator<< (ostream &os, Well &well);
  
  friend class  P0Eqn;
  friend class CYPEqn;
  friend class CPPEqn;
  friend class  P2Eqn;

private:
  bool debug;

  int wLens; 
  int    *wLocIJK[3];   // multiblock penetration 
  double *wLocXYZ[3];    // xyz location 
  double *wTrans, wStr;  // rate for fixed rate, and pressure for fixed P 
  WConType wConType;

  double calcQ();        // calc q for BHP controlled well
  double calcPw();       // calc BHP for rate controlled well 
  void convertUnit();    // convert Field unit to Matric unit

  bool isProdWell;
  bool isInjcWell;
  double angle[2];
  double radius;
  int npts;

  // calcQ(), calcPw(), update() are not used now, and not implemented
};

#endif
