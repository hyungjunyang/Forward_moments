//////////////////////////////////////////////////
//                                              //
//                 Liyong Li                    //
//      Reservoir Simulation Research Team      //
//        Chevron Petroleum Technology Co.      //
//          1300 Beach Blvd., Rm. 3166          //
//            La Habra, CA 90631-6374           //
//             Phone:(562) 694-7366             //
//              Fax:(562) 694-7565              //
//            Email liyl@chevron.com            //
//                                              //
//                 Version 1.                   //
//                Apr. 07, 1999                 //
//////////////////////////////////////////////////

//             --       MCTrvl.h          --
//             --    class definition    --

#ifndef MCTrvlH
#define MCTrvlH

#include <iostream>
#include <math.h>
#include <assert.h>
#include <fstream>
#include <iomanip>

class MCTrvl{
public:
   MCTrvl();
  ~MCTrvl();
   int getLength() {return NLn;}
   double * getQtAvg() {return qt_avg;}
   double * getQt2Avg() {return qt2;}
   double * getTrAvg() {return tr_avg;}
   double * getTrVar() {return tr_var;}
   double * getRho() {return r_tau1_tau2;}
   
// OVERLOADED OPERATORS:
// MCTrvl& operator= (const MCTrvl &rhs);
   void printTrvlTime();
   private:    
   int NLn, Nr;
   double r_given;
   double *theta, *thetaw;
   double *tr_rand, *v1_rand, *v2_rand;
   double *tr_avg, *tr_var;
   double *qt_avg;
   double *r_tau1_tau2;
   double *qt2;

//   double *th1, *th2;
//   double *q1_avg, *q1_var;
//   double *q2_avg, *q2_var;

   void initialize();
   void readdata();
   void adjust_flux();
   void calcVeloMoments();
   void calcLogTrvlMoments();
   void calcTrvlMoments();
   void calcStatistics(int&, double *, double&, double & );
   void test_Trvl_Distribution();
   void normal2lognormal(double&  tauAvg, double&  tauVar, 
                         double& ltauAvg, double& ltauVar); 
   void normal2lognormal(double& tauAvg1, double& tauVar1,
                         double& tauAvg2, double& tauVar2, double& ctau1tau2,
            double& ltauAvg1, double& ltauVar1, double& rho12);
};

#endif
