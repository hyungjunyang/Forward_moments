#ifndef _GaussCovariance_H
#define _GaussCovariance_H
#include <math.h> 
#include "Common.h"
#include "Covariance.h"
#include <iostream>

class GaussCovariance : public Covariance{
  public:
   GaussCovariance(double var, double cor1, double cor2, double cor3, double azimuth);
   virtual ~GaussCovariance();
   virtual double getCov(double x1, double y1, double z1, 
                         double x2, double y2, double z2);
  private:
};

#endif
