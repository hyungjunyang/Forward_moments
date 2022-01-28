#ifndef MC_TO_SME_H
#define MC_TO_SME_H

#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<math.h>

class MC_TO_SME{

public:
  MC_TO_SME(int);
 ~MC_TO_SME();

void initialize();
void read();
void print();
void transform1(int& len, double* tau, double *var, double* lntau, double* lnvar);
void transform2(int& len, double* tau, double *var, double* tau1tau2, double* rho);

double* getAvg()    {return avg;}
double* getVar()    {return var;}
double* getCov()    {return cov;}
double* getLogAvg() {return avg_log;}
double* getLogVar() {return var_log;}
double* getLogCov() {return cov_log;}

private:
  int length;
  double *avg    , *var    , *cov;
  double *avg_log, *var_log, *cov_log;
};

#endif
