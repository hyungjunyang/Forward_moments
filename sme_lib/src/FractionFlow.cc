#include "FractionFlow.h"

FractionFlow::FractionFlow(double Swc, double Sor, double viscosR)
:Swc_(Swc), Sor_(Sor), viscosR_(viscosR) 
{        
   deltaS = 1.0 - Swc_ - Sor_;
   sStar = Swc_ + fabs( deltaS ) / sqrt(viscosR_ + 1);
   fpwStar = fpw(sStar);
}

FractionFlow::~FractionFlow() {
}
/////////////////////////////////////////////////////////
//                                                     //
//        Function: Saturation in Lagrangian Space     //
//                  or the deterministic solution      //
//                                                     //
//        Corey-Type relative permeability is used     //
//                                                     //
//                 mS*S                                //
//        fw=----------------------                    //
//            m*S*S+(1-S)*(1-S)                        //
//                                                     //
//            MUn                                      //
//        m = ----                                     //
//            MUw                                      //
//                                                     //
//                2mS*(1-S)                            //
//        f'w=----------------------; expressas fpw    //
//              [m*S*S+(1-S)*(1-S)]^2                  //
//                                                     //
//                     1                               //
//        S* = -------------------                     //
//                 SQRT(1 + m)                         //
//                                                     //
/////////////////////////////////////////////////////////
// --- function to calculate fw and fpw at a given sw ---
// --- (1) Corey Type relative Permeability ---
double FractionFlow::fw(double sw) {
  return viscosR_ * (sw - Swc_) * (sw - Swc_)
      /( viscosR_ * (sw - Swc_) * (sw - Swc_) + (1. - sw - Sor_ )*(1. - sw - Sor_) );
}

//      d fw
//      ----
//      d sw
double FractionFlow::fpw(double sw)   {
   return (2. * viscosR_ * deltaS * (sw - Swc_) * (1. - sw - Sor_) )
    /pow(viscosR_ * (sw - Swc_) * (sw - Swc_) + (1. - sw - Sor_)*(1. - sw - Sor_), 2);
}
/*
////////  Tracer Case ////////////////////////////////////
//                                                      //
//                                                      //
//        fw = mS         and M = 1                     //
//                                                      //
double  FractionFlow::fw(double sw) {
   return sw - Swc_; 
}
double FractionFlow::fpw(double sw) {
   return 1.;
}
//                                                      //
//////////////////////////////////////////////////////////
*/
