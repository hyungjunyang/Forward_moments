/*
 * meas -- measurement
 * uncd -- unconditional
 * cond --   conditional
 * 
 */

#ifndef _PRESSMEASURE_H
#define _PRESSMEASURE_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <math.h>
#include "Grid.h"

class PressMeasure {
public:
  PressMeasure(int num, Grid *, bool);
 ~PressMeasure();
 
  //set fuctions
  void setPdata(int m, int i, int j, int k, double p);
  void setPcond(int m, double p, double p_std);
  void setPuncd(int m, double p, double p_std);
  void setPcond(int m, double p);
  void setPuncd(int m, double p);

  // get functions 
  int  getLength()      { return num_p_meas;}
  int  getIpmeas(int i) { return i_pmeas[i]; }
  int  getJpmeas(int i) { return j_pmeas[i]; }
  int  getKpmeas(int i) { return k_pmeas[i]; }
  int  getIJK(int i) { return ijk_pmeas[i]; }

  int* getIpmeas()      { return i_pmeas; }
  int* getJpmeas()      { return j_pmeas; }
  int* getKpmeas()      { return k_pmeas; }
  int* getIJKpmeas()    { return ijk_pmeas; }

  double  getPmeas(int i)  { return p_meas[i];}
  double  getPpred(int i)  { return p_avg_cond[i];}
  double  getWpmeas(int i) { return w_pmeas[i];}
  double* getPmeas()  { return p_meas;}
  double* getPpred()  { return p_avg_cond;}
  double* getPpstd()  { return p_std_cond;}
  double* getWpmeas() { return w_pmeas;}

  double getAbsDiff(double&);
  double getAbsDiff_D(double &);

  double getSqdDiff();

  double getAbsDiff(double&, double&); 
  double getSqdDiff(double&);
  
  double* getPStdCond() {return p_std_cond;}
  double* getPStdUnCd() {return p_std_uncd;}
  double* getPAvgCond() {return p_avg_cond;}
  double* getPAvgUnCd() {return p_avg_uncd;}  
  
  // output
  void output();
  void output(char *);
  void display();
  
private:
  bool debug;
  Grid *gridPtr;
  int num_p_meas;
  int    *i_pmeas, *j_pmeas, *k_pmeas, *ijk_pmeas;
  double *p_meas;
  double *w_pmeas;

  double *p_avg_cond;
  double *p_avg_uncd;
  double *p_std_cond;
  double *p_std_uncd;
  
  void initilize();
  double abs_avg_diff;
  double abs_std_diff;
};
#endif
