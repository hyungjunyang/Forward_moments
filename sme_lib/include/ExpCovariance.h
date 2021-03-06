#ifndef _ExpCovariance_H
#define _ExpCovariance_H
#include <math.h> 
#include "Common.h"
#include "Covariance.h"
#include <iostream>

class ExpCovariance:public Covariance{
  public:
   ExpCovariance(double var, double cor1, double cor2, double cor3, double azimuth);
   virtual ~ExpCovariance();
   virtual double getCov(double x1, double y1, double z1, 
                         double x2, double y2, double z2);
  private:
};

#endif
