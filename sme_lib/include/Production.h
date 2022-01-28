#ifndef _PRODUCTION_H_
#define _PRODUCTION_H_

#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>
#include <time.h>

#include "Common.h"
#include "Control.h"
#include "Saturation.h"
#include "BLSolution.h"
#include "MCTrvl.h"
#include "Well.h"
#include "SwCondiInverse.h"

double pdfGauss(double& x, double& avg, double& var);
double pdfGauss2(double& x1, double& avg1, double& var1,
		 double& x2, double& avg2, double& var2, double& rho
		);
double pdfGauss2(double& x1, double& x2, double& rho);

class Production {
  public:
  Production(Unit u, BLSolution* blSoln, int, double, bool debug_);
  Production(Unit u, Domain* domain, BLSolution* blSoln, int, double, bool debug_);
  Production(Unit u, Domain* domain, Flow* flow, BLSolution* blSoln, int, double, bool debug_);
 ~Production();
 
  void calcWaterCutMomentFast(int, double*, double*, double*);
  void calcWaterCutMomentFast(int, double*, double*, double*, double*);
  void calcWaterCutMoments(int, double*, double*, double*);
  void calcWaterCutMoments(int, double*, double*, double*, double*);
  void calcWaterCutMoments(int, double*, double*, double*, double*, double*); 
  void calcWCumMoments(int, double*, double*, double*, double*, double*);
                                 
  void calcWaterCut(ifstream &is, char *Dir);
  void convertUnit();

  // --- Functions for ProdCondiInverse (Pipatl)
  double CalcTotalProdRateMean(int wInd);
  double CalcTotalProdRateCov(int wInd1, int wInd2);
  double CalcWaterProdRateMean(int wInd, double wTime);
  double CalcWaterProdRateCov(int wInd1, double wTime1, double qwMn1, int wInd2, double wTime2, double qwMn2);
  double CalcTotal_WaterProdRateCov(int wTInd, int wwInd, double wwTime);
  double CalcCrossCovQwY(int wInd, double wTime, int yInd);

  // --- Functions for Water Production (Pipatl)
  void ProdWellIs(int wInd);
  void AddAllProdWells();
  void SetParticlesAtProdWells();
  void LaunchParticles();
  void CalcCovYLnTau(int wInd);
  void ClearCovYLnTau();

private:
  bool debug;
  Unit unit;

  Domain* domainPtr;
  BLSolution* blSolnPtr;
  ParticleTrack* ptclTrackPtr;
  Point* pt_start;
  Grid* gridPtr;
  Flow* flowPtr;

  int     num_steps;
  double  t_maximum;
  double  dt;
  double *time_step;
  double *water_production;

  // launch particles for water producers
  static const int int_section_per_side = 6;
  static const int int_quad_point = 3;  
  int nWaterProd;
  int nParticles;
  map<int,int> prodWellMap;
  struct WellParticles{
    int start;
    bool side[6];
  };
  WellParticles* wParticlesPtr;
  Point* pointPtr;
  ParticleTrack* bTrackPtr;

  // Calculate water production rate
  double** covYLnTau;
  double** covVY;

  // water cut and oil cut
  double *w_cut_avg, *w_cut_var;
  double *o_cut_avg;
  double *w_cut_var2; // q^2(t,t') 
  double *w_cut_var_1, *w_cut_var2_1;

  // cumulative quantities
  double *w_cum_avg, *w_cum_var;
  double *o_cum_avg, *o_cum_var;

  // numerical integration (Gaussian Quadrature)
  int num_intg;
  double *tau,  *tau_weight,  *fw;
  double *tau1, *tau_weight1, *fw1;
  double *xx, *xw;
  double *gx, *gw, *xwxw;
  double multiple;
  // --- Private Functions --- 
  void initialize();
  void time_discretize();
};

#endif
