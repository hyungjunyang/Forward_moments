#include "Covariance.h"

Covariance::Covariance(double var, double cor1, double cor2, double cor3, double azimuth_)
: variance (var), cor_length1(cor1), cor_length2(cor2), cor_length3(cor3), azimuth(azimuth_)
{
}

Covariance::~Covariance(){
}

 
