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

//          --       Particle.h       --
//          --    class definition    --

#ifndef PARTICLEH
#define PARTICLEH

#include "Domain.h"
#include <vector>


double pdfLognormal(double& x, double& avg, double& var);

class Particle{
	
public:
   Particle(Point p, double time_, Domain&);
   virtual ~Particle();

   void setBlock(Block*); 
   Block* getBlock();	
   void move();
   void moveBackward();//By Pipat
   void reorderData();//By Pipat

   void initialization();
   void CoordTrans();
   void CoordTrans2(double *Vxi_Vxi, double *Vxi_Eta, 
		    double *Eta_Eta);
   void calcTauVar();
   
   void calcTauAvg();
   void calcTauAvg0();
   void calcTauAvg2();
   
   void calcTauAvg2Trans(double*, double*);
   
   void calcTauAvgPdf();
   void calcEtaDU(double *etadu);
 
   double getTrvlTimeAvg(int i)  {return tau_Avg0Trans[i]; }
   //double getTrvlTimeAvg(int i)  {return tau_Avg2Trans[i]; } //2nd order
   double getTrvlTimeVar(int i)  {return tau_VariTrans[i]; }
   
   double getLastTravelTime() {return getTrvlTimeAvg(length-1);}
   double getLastTravelTVar() {return getTrvlTimeVar(length-1);}
   double getLastTravelVBar() {return v_avg[0];}
   int    getLength()         {return trajectory.size();}
   double getDs(int i)   {return ds[i];}
   double getVAvg(int i) {return v_avg[i];}
   double getDvDeta(int i) {return dv_deta[i];}
   Point& getTrajectory(int i) {return trajectory[i];} 
   Point& getTrajectory()      {return getTrajectory(getLength()-1);}
   Point  getParticle_Point(){return particle_point;}//By Pipat for testing backtracking only
   vector<Point> getVecTrajectory(){return trajectory;}//By Pipat for testing backtracking only
   vector<double> getVecTimeFlight(){return time_flight;}//By Pipat for testing backtracking only
   int    getNode(int i) {return node[i];}
   int    getNode_i(int i) {return node_i[i];}
   int    getNode_j(int i) {return node_j[i];}
   double getCosth(int i) {return costh[i];}
   double getSinth(int i) {return sinth[i];}
   double getX_x1(int i)  {return x_x1[i];}
   double getX2_x(int i)  {return x2_x[i];}
   double getY_y1(int i)  {return y_y1[i];}
   double getY2_y(int i)  {return y2_y[i];}

   // ===== Cartesian System =====   
   void try_Carte();
   void calcTauAvg0Carte(double*, double*);
   void calcEtaEtaCarte(double*, double*, double*, double*, double*,
		        double*, double*, double*, double*, double*, double*);
   void calcTauAvg2Carte(double*, double*, double*, double*,
		         double*, double*, double*, double*);
   void calcTauVariCarte(double*, double*, double*, double*, double*, double*);
                         
   void PrintTauCarte(ostream & os);
   // ===== end of Cartesian System =====
   
   void PrintParticle(ostream & os);
   void PrintPathline(ostream & os);
   void PrintTrajectory(ostream & os);
   void PrintVelocity(ostream & os);
   void PrintTrajectory();

private:
   Point  particle_point;
   Block* particle_block;
   Domain& domain_;
   double particle_time_;
   vector<Point> trajectory;
   vector<double> time_flight;
   
   int length;
   int *node, *node_i, *node_j;
   double *costh, *sinth;
   double *ds; 
   double *dv_deta;
   double *x_x1, *y_y1, *x2_x, *y2_y;

   double *v_avg, *v_var;
   
   // ===== Original System =====
   double *tau_Avg0;
   double *tau_Avg2;
   double *tau_Vari;
  
   // ===== Transformed System =====
   double *tau_Avg0Trans;
   double *tau_Avg2Trans;
   double *tau_VariTrans;

   // ===== Cartesian System =====   
   double *tau_Avg0Carte;
   double *tau_Avg2Carte;
   double *tau_VariCarte;

   // ===== end of Cartesian System =====

   double *tau_AvgPdf;
   bool try_Cartesian;
   
};

#endif
