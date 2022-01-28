//////////////////////////////////////////////////
//                                              //
//                 Liyong Li                    //
//      Reservoir Simulation Research Team      //
//      ChevronTexaco Energy Technology Co.     //
//        Bollinger Canyon Rd, Rm. D2102        //
//            San Ramon, CA 94583               //
//             Phone:(925) 842-6359             //
//              Fax :(925) 842-6283             //
//         Email liyl@chevronTexaco.com         //
//                                              //
//                 Version 1.                   //
//               July. 14, 2004                 //
//////////////////////////////////////////////////
//                                              //
//      Based on Dagan's Paper 1985             //
//                                              //
//      p_avg_UnCd                              //
//      p_std_UnCd                              //
//      p_avg_Cond                              //
//      p_std_Cond                              //
//      p_avg_Filt                              //
//      p_std_Filt                              //
//                                              //
//////////////////////////////////////////////////
//                                              //
//          --    PCondiInverse.h          --   //
//          --    class definition    --        // 
//////////////////////////////////////////////////


#ifndef _PCondiInverse_H
#define _PCondiInverse_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <iomanip>
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
#include "KrigWeigh.h"

class PCondiInverse {
public:
 // PCondiInverse(Grid*, Perm*, const P0Eqn &P0ee);
  PCondiInverse(Grid*, Perm*, P0Eqn *P0ee, CYPEqn *CYPe, CPPEqn *CPPe, bool);

 ~PCondiInverse();
 
  void solve(int nWells, const Well *wellArr, const Control &c, 
             Solver &s, int time_step );

  void solveIter(int nWells, const Well *wellArr, const Control &c, 
                 Solver &s, int time_step );
	  

  //output()
  void output_p(char*, double*, double*);
  void output_p_Meas(char *file);

private:
  bool debug;

  // object Pointers
  Grid   *gridPtr;
  Perm   *permPtr;
  P0Eqn  *P0e;
  CYPEqn *CYPe;
  CPPEqn *CPPe;
  const double *p0;
  
  // krigging weight
  KrigWeigh      *krigPtr;
  void delKrigWeighObj();
  void addKrigWeighObj();

  // pressure measurements
  int num_p_Meas, meas_err_opt;
  int *i_p_Meas, *j_p_Meas, *k_p_Meas, *ijk_p_Meas;
  double *p_avg_Meas, *p_std_Meas;
  double *meas_err_mat;

  // pressure UnConditional and Condtional
  double *p_avg_UnCd, *p_std_UnCd;
  double *p_avg_Cond, *p_std_Cond;
  double *p_avg_Filt, *p_std_Filt;
 
  // permeability UnConditional and Condtional
  double *Y_avg_UnCd, *Y_std_UnCd;
  double *Y_avg_Cond, *Y_std_Cond;
  
  // Iteration Related Fields
  int max_iter;
  double eps1, relax;
  double abs_avg_diff;
  double abs_std_diff;
  double *error_old;
  double *error_crt;
  double *error_new;
  void setIterFields();
  void delIterFields();
  double getDiff(double *, double *);
  void calcMu(double*, double*, double* mu);
   
  void condCPPbyGaussDist();
  void condCPPbySolveEqn(int nWells, const Well *w, 
		         const Control &c, Solver &s);
  
  // private functions
  void readData();
  void addMem_of_readData();
  void delMem_of_readData();
  void addArrays();
  void delArrays();
  void setUnCondMoments();
};
#endif
