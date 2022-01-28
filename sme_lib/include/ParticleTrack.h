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

//          --     ParticleTrack.h     --
//          --    class definition    --


#ifndef PARTICLETRACKH
#define PARTICLETRACKH

#include "Domain.h"
#include "Particle.h"
#include "WellArc.h"
#include "Point.h"
#include "Common.h"

void transform1(int& len, double* tau, double *var, double* lntau, double* lnvar);
void transform2(int& len, double* tau, double *var, double* tau1tau2, double* rho);

class ParticleTrack{ 

public:
   ParticleTrack(Domain* domainPtr, int nstrl, Point* pt,
                 bool print_option = false);
   
   ParticleTrack(Domain* domainPtr, bool debug, bool print_option = false);
   ~ParticleTrack();

   void initialization();

   void goTrack();
   void goTrackBackward();//By Pipat
   void calcTrvlTimeMoments();
   void calcTrvlTimeMomentsNoWell();//By Pipat
   
   void calcTrvlTimeCorrelation();
   void calcCqqAtWell(double, double*);

   void calcMoments(Particle*, Particle*, double *Vxi_Vxi, double *Vxi_Eta,
                    double *Eta_Eta, double *, double*);
   void calcMoments(Particle*, Particle*, double *Vxi_Veta, double *Vxi_Eta);
                    
   int     getNumStrl()         {return nstrl_;}
   double* getLogTrvlTimeAvg()  {return trvl_avg_log;}
   double* getLogTrvlTimeVar()  {return trvl_var_log;}
   double* getLogTrvlTimeRho()  {return trvl_rho_log;}   
   double* getVbar()            {return vbar;}
   double* getCqq()             {return Cqq;}

   Particle* getParticle(int i) {return particles[i];}
   
   void printTravelTime();
   void printLogTravelTime();
private:
   bool debug;

   Unit unit;
   Domain* domainPtr;
   WellArc* prodWellArcPtr;
   Block * prodBlkPtr;

   // Cv1v1
   double *Cv1v1, *Cv2v2, *Cv1v2;

   int nstrl_;
   Particle**  particles;

   // travel time mean, variance, and correlation;
   double *trvl_avg    , *trvl_var    , *trvl_rho;
   double *trvl_avg_log, *trvl_var_log, *trvl_rho_log;
   double *vbar;
   double *Cqq;

   // print tag
   bool print_status;
};

#endif
