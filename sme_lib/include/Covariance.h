#ifndef _Covariance_H
#define _Covariance_H
#include <math.h> 
#include "Common.h"

class Covariance{
  public:
  Covariance(double var, double cor1, double cor2, double cor3, double azimuth_);
  ~Covariance();
   virtual double getCov(double x1, double y1, double z1, 
                         double x2, double y2, double z2) = 0;
  protected:
   double variance, cor_length1, cor_length2, cor_length3, azimuth;
};

#endif
