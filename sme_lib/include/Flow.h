#ifndef _FLOW_H_ 
#define _FLOW_H_

#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <time.h>

#include "Common.h"
#include "Domain.h"
#include "Control.h"
#include "Grid.h"
#include "Region.h"
#include "Perm.h"
#include "Fluid.h"

#include "Eqn.h"
#include "P0Eqn.h"
#include "CYPEqn.h"
#include "CPPEqn.h"
#include "CPPDEqn.h"
#include "P2Eqn.h"
#include "Well.h"
#include "Solver.h"
#include "FullLUSolver.h"
#include "LapackBandSolver.h"
#include "LapackFullSolver.h"
#include "LapackSymBandSolver.h"
#include "Velocity.h"

using namespace std;

class Flow {
public:
  Flow(char *Dir, Domain *, bool); 
 ~Flow();
  Velocity* getVeloPtr() {return veloPtr;}
  Perm* getPermPtr() {return permPtr;} // Pipatl
  Eqn* getP0e() {return P0e;} // Pipatl
  Eqn* getCYPe() {return CYPe;} // Pipatl

  void solve(char *Dir, ifstream &is);
  void writeToFile(char *Dir);
  void writeToDomain(Domain * domainPtr);
  void deleteVeloPtr(){delete veloPtr;} // (Pipatl) This is to be used before calling solve() repeatedly

  //Flow(ifstream &is, char *Dir);
  //Grid* getGrid()        {return gridPtr;}
  //int getNumOfWells()    {return nWells;}
  //Well* getWellArr()     {return wellArr;}

  //void testBackTracking(){SwCondiInverse* swcondiPtr = new SwCondiInverse(gridPtr,domainPtr);delete swcondiPtr;}//By Pipat test backtracking

private:
  bool debug; 
  int pressure_condition;

  void initialize(char *Dir);

  Unit     unit;
 
  Domain* domainPtr;

  //void readUnit  (ifstream &is );
  Grid *gridPtr;
  //void addGridObj(ifstream &is, char * Dir);

  Region *regnPtr;
  //void addRegnObj(ifstream &is, Unit unit);

  Perm *permPtr;
  void addPermObj( Unit unit, char *Dir);
  
  Fluid *fluidPtr;
  //void addFluidObj(ifstream &is, Unit unit );
  
  int nWells;
  Well    *wellArr;
  //void  addWellArr(ifstream &is );

  int nDt;
  Control *contPtr;
  //void  addContObj(ifstream &is );
  
  Eqn *P0e, *CYPe, *CPPe, *P2e;
  void  addPresObj();
  
  Solver *ss;
  void addSolverObj();
  void delSolverObj();

  Velocity *veloPtr;
  void   addVeloObj( char *Dir );
  

};

#endif
