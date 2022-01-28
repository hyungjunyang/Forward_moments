/*
 * File: Saturation.h
 * --------------------------------------------------------------------
 * Saturation moments calculation
 *   1) input travel time mean and variance
 *   2) calculate saturation mean and variance
 *
 * --------------------------------------------------------------------
 **Streamline Data***************************************************** 
 *   int nLn -- number of Streamlines,  
 *   int nNd -- number of point at each streamline, 
 *   int nPoint = nLn * nNd
 *   double *tauMn, *tauVari -- travel time mean and variance 
 *   double *x1C, *etaMnC    -- streamline coord in old Cartesian grid
 *                                    needed on visualization
 * 
 **Saturation Data*****************************************************
 *
 *   double sStar, fpwStar -- front saturation and velocity
 *   double visR,          -- viscosity ratio (mu_o/mu_w)
 *   double Swc, Sor,
 *   deltaS                -- end point saturation, deltaS=1-Swc-Sor 
 *
 *   int nS          -- Number of points to discrete  saturation profile  
 *   double *tau     -- the array of tau to calculate saturation profile
 *   double *s       -- the discreted saturation profile(tau vs. sat) 
 *   
 *   double t        -- require solution at time t
 *   double *satDet  -- deterministic saturation solution at time t
 *   double *satAvg  -- saturation mean at time t
 *   double *satVar  -- saturation variance at time t
 *
 * --------------------------------------------------------------------
 * Created    08/03/00          hcao     
 * Modified   09/04/04          liyl
 */

#ifndef _SATURATION_H
#define _SATURATION_H

#include <fstream>
#include <math.h>
#include <time.h>
#include "Common.h"
#include "ParticleTrack.h"
//#include "pdf.cc"
//#include "quadrature.cc"

//added by liyl<<
double pdfGauss(double& x, double& avg, double& var);
double pdfGauss2(double& x1, double& avg1, double& var1, double& x2, double& avg2, double& var2, double& rho);
double pdfGauss2(double& x1, double& x2, double& rho);
//added by liyl>>

void GaussianQuadrature(int qPoint, double* &rloc, double* &weight);// Pipiatl or quadrature.cc

class Saturation {

public:
  Saturation(ParticleTrack *ptclTrackPtr_, 
   double visR_, double Swc_, double Sor_, int nS_ );
  ~Saturation();

  void solve(double tt);                                // calculate saturation moments
  void output(int i);									// output

  void calcSatCovEndPoints();	// calc sat cov among end points of all streamlines (Pipatl)
  double* getSatCovEP(){return satCovEP;}	//Pipatl
  double getSatAvg(int nStreamline,int nPoint){return satAvg[(nStreamline-1)*nNd+nPoint-1];}//Pipatl
  double getSatVar(int nStreamline,int nPoint){return satVar[(nStreamline-1)*nNd+nPoint-1];}//Pipatl


private:
  friend class SwCondiInverse;
  bool debug_;

  ParticleTrack *ptclTrackPtr;

  // --- Streamline Data ----------------------
  int nLn, nNd, nPoint;
  double *tauMn, *tauVar;
  double *tauRhoLog; //Need to calculate Css (Pipatl)
  double *x1C, *etaMnC;

  // --- Saturation data ---------------------
 
  double sStar, fpwStar, visR, Swc, Sor, deltaS; 
  double *satDet, *satAvg, *satVar;
  double *satCovEP; //Covariance of saturation at end points of all streamlines(Pipatl)
  double  t;
   
  int nS;
  double *tau, *s;
  
  void initialize();                                     // initialize
  void calcSatDet(int n, const double *ta, double *ss);  // calc deter sat
  void calcSatAvgVar();                                  // calc sat mean,vari
  void calcSatAvgVar2();								 // calc sat mean,vari by dlntau(Pipatl)
  void calcSatAvgVarQuad();								 // calc sat mean, vari by quadrature intregation(Pipatl)
  double calcSatCov(double satAvg1, double satAvg2, double tauMn1, double tauMn2, double tauVar1, double tauVar2, double tauRhoLog12);
  // calc sat covariance by integration over bivariate pdf for point 1 and 2 (Pipatl)
  double calcSatCov2(double satAvg1, double satAvg2, double tauMn1, double tauMn2, double tauVar1, double tauVar2, double tauRhoLog12);
  // calc sat covariance by integration over bivariate pdf for point 1 and 2 using dlntau(Pipatl)
  double calcSatCovQuad(double satAvg1, double satAvg2, double tauMn1, double tauMn2, double tauVar1, double tauVar2, double rhoLntau12);
  // calc sat covariance by quadrature integration over bivariate pdf for point 1 and 2 (Pipatl)

  double fw(double sw);                                  // fw(sw)
  double fpw(double sw);                                 // Dfw_Dsw(sw)
  double root(double a, double eps);                     // root of fpw-a=0
  double pdfForTau(double ta, double tMn, double tVari); // tau's pdf
};

#endif

