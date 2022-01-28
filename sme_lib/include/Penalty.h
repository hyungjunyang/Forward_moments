/* 
num_mpts : the number of the master points.
cpq      : the C matrix in the gradient equation 
ppq      : the d vector in the gradient equation
xpq      : the dY, or the perturbations of Y
one_neg  : -[ I ], the negative identity matrix
one_pos  :  [ I ], the positive identity matrix
bpq_min  : defined as Y - Y_min
bpq_max  : defined as Y_max - Y
*/

#ifndef _PENALTY_H
#define _PENALTY_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <math.h>

#include "Perm.h"
#include "MasterPoint.h"
#include "PressMeasure.h"
#include "P0_Sensitivity.h"

class Penalty {
public:
   Penalty(Perm *perm, MasterPoint* mpts, PressMeasure *prsmPtr, 
           P0_Sensitivity *sensPtr);
  ~Penalty();   
  
   void setP0data();
   void setY0data();
   void setYConddata();

   double getAvrgDiff() { return avrg_diff; }
   double getFobjNorm() { return fobj_norm; }
   double getFobjIter() { return fobj_iter; }
   double getFobjInit() { return fobj_init; }
   
   void buildMatrix();
   void minimize();
   void calcGrad();
   void calcNumActRst();
   void projectGrad();
   
   void initObjValue();
   void calcObjValue();
   void iterObjValue();

   void output();
private:
   Perm         *permPtr; 
   int num_mpts;

   MasterPoint  *mptsPtr;

   int num_p_meas;
   PressMeasure *prsmPtr;

   P0_Sensitivity *sensPtr;
   
   double fobj_init;
   double fobj_iter;
   double fobj_norm;
   double avrg_diff, avrg_diff2;

   double *cpq, *ppq, *xpq;
   double *one_neg, *one_pos;
   double *bpq_min, *bpq_max; 

   double *p_pred_meas;
   double *p_pred_meas_wpms;

   double *Y_avg,      *Y_var;
   double *Y_avg_cond, *Y_var_cond;

   int num_active;
   int    *index_active;
   double *grad;
   double *u, *uprim, *hh_t1;

   void initilize();
};

#endif


