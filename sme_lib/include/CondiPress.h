#ifndef _CONDI_PRESS_H
#define _CONDI_PRESS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <stdlib.h>
#include <math.h>

#include "Grid.h"
#include "Perm.h"
#include "Eqn.h"
#include "P0Eqn.h"
#include "CYPEqn.h"
#include "CPPEqn.h"
#include "Control.h"
#include "Solver.h"
#include "Well.h"

#include "MasterPoint.h"
#include "PressMeasure.h"
#include "KrigWeigh.h"
#include "P0_Sensitivity.h"
#include "PVar_Sensitivity.h"

#include "Penalty.h"

using namespace std;

class CondiPress {
public:
 // CondiPress(Grid*, Perm*, const P0Eqn &P0ee);
  CondiPress(Grid*, Perm*, P0Eqn *P0ee, CYPEqn *CYPe, CPPEqn *CPPe, bool);

 ~CondiPress();
  void solve(int nWells, const Well *w, const Control &c, 
             Solver &s, int );
  
  void solve_pstd(int nWells, const Well *w, const Control &c, 
                  Solver &s, int );
  void solve_p0(int nWells, const Well *w, const Control &c, 
                  Solver &s, int );
  void calcSensitivity(int nWells, const Well *w, 
                       const Control &c, Solver &s, double*, double*);

  //output()
  void output(char*, double*, double*);
  void output();
  void output_p0();
  void outputCondiPerm();

  // try()
  void try_Dagan(int nWells, const Well *wellArr, const Control &c, 
                           Solver &s, int time_step );
private:
  bool debug;

  Grid         *gridPtr;
  Perm         *permPtr;
  //const P0Eqn  *P0e; 
  P0Eqn  *P0e;
  CYPEqn *CYPe;
  CPPEqn *CPPe;
  const double *p0;
  
  // master point
  MasterPoint *mptsPtr;
  void addMasterPointObj(long seed, int num_i, int num_j, int num_k);

  // krigging weight
  KrigWeigh      *krigPtr;
  void addKrigWeighObj();

  // pressure measurements
  PressMeasure *prsmPtr, *pstdPtr;
  void addPressMeasureObj(int num);

  // pressure sensitivity
  P0_Sensitivity *p0snPtr;
  void addP0SensitivityObj();
 
  // pressure variance sensitivity
  PVar_Sensitivity *pVsnPtr;
  void addPStdSensitivityObj(int );
  void delPStdSensitivityObj();  
	  
  // Penalty function
  Penalty        *pnltPtr;
  void addPenaltyObj();  

  // penalty and minimization
  int num_p_meas, num_mpts;
  int num_pstd_meas;
  int max_iter, min_iter, num_perm_meas, num_step;
  int ind_plot, ind_fobj, ind_dbug, ldbg;
  double eps1, eps3, eps5, relax, dconve, multplier;
  
  double *cpq;
  double *apq;
  double *bpq;
  double *ppq;
  
  double *cpq_var;
  double *bpq_var;
  double *ppq_var;

  double *Y_avg_UnCd, *Y_std_UnCd;
  double *Y_avg_Cond, *Y_std_Cond;
  double *Y_avg_Mast, *Y_std_Mast;
  
  double *p_avg_UnCd, *p_std_UnCd;
  double *p_avg_Cond, *p_std_Cond;
  double *p_avg_Meas, *p_std_Meas;
  double *p_avg_Old,  *p_std_Old;
  double *p_avg_New,  *p_std_New;
  
  double *dY_avg, *dY_std;
  double *dpAvg_dY, *dpStd_dY;

  double sumdiv, sumwall;
  double fobj0, fobj_ini;
  double sumdiv_var;
  double fobj0_var, fobj_ini_var;
 
  int num_iter_control;

  // private functions
  void readData();
  void initialization();
  
  void set_P0_in_prsmPtr();
  void set_P0_in_prsmPtr(double*, double*);
  void set_PStd_in_pstdPtr();
  
  void setYCondAtMstPts();
  void setYUnCdAtMstPts();
  
  void addLocalArrays();
  void delLocalArrays();
  void setUnCondMoments();

  void condCPPbyGaussDist();
  void condCPPbySolveEqn(int nWells, const Well *w, 
		         const Control &c, Solver &s);
  void krigAtMasterPoints();

  void dYWithrelax();
  void calcPenaltyF();
  void calcPenaltyFD();
  void calcPenaltyFVar();
  void minimization();
  void minimizationVar();

  // Temporal Try Area
  /*
   int num_temp;
   int *ii_tmp, *jj_tmp, *kk_tmp;
   */
};
#endif
