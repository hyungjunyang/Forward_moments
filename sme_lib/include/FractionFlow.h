#ifndef _FRACTIONFLOW_H
#define _FRACTIONFLOW_H
#include <math.h>

class FractionFlow {
public:
   FractionFlow(double, double , double ); 
  ~FractionFlow();
  double getSStar()   { return   sStar; }
  double getFpwStar() { return fpwStar; }
  double fpw(double sw);   // Dfw_Dsw(sw)
  double fw(double sw);                                  // fw(sw)

private:
  double sStar, fpwStar;
  double deltaS;
  double viscosR_;
  double Swc_, Sor_;
};

#endif
