#ifndef _KRIGWEIGH_H
#define _KRIGWEIGH_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <math.h>

#include "Grid.h"
#include "Perm.h"

class KrigWeigh {
public:
   KrigWeigh(Grid*, Perm*, int num);
  ~KrigWeigh();

   void initialization();
   int getIndex(int& i, int& j) { return i + j * num_meas;}
   int getLength() {return length; }
   void calCovInverse(int meas_len, double* cov22, double* cov22inv);
   void calWeigh(int* im, int* jm, int* km);
   void condGaussDist(int cova_len, int meas_len, int* node, 
                      double *cova, double*);

   double* getWeigh() {return weigh;} 
   double  getWeigh(int& i, int& j) {return weigh[getIndex( i, j) ];}
   void display(); 
   void output();
   
private:
   bool debug;
   Grid* gridPtr;
   Perm* permPtr;

   int num_meas, length;
   double* weigh;

   void setCovariance22(int* im, int* jm, int* km, double* cov22);
   void setCovariance12(int* im, int* jm, int* km, double* cov12);
   void calCovInverse(double* cov22, double* cov22inv);
   void linv(int& n, double* a, double* b);
   int  cholsky(int& n, double* a, double* t);
};

#endif
