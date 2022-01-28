/*
 * File: BLSolution.cc
 * ----------------------------------
 * Implementation for BLSolution class
 */
#include "BLSolution.h"

using namespace std;

// --- Constructors ---
BLSolution::BLSolution(int nPts, double Swc, double Sor, double viscosR)
:nPts_(nPts), Swc_(Swc), Sor_(Sor), viscosR_(viscosR) {
  fractionFlow = new FractionFlow(Swc_, Sor_, viscosR_);
  initialize();
  cal_BL_solution();
}

BLSolution::~BLSolution() {
  delete[] bl_1d;
  delete[] tau;
  delete[] bl_slope;
  delete fractionFlow;
}

// --- initialize BLSolution data ----------------
void BLSolution::initialize() {
  eps = 1.0E-12;
  sStar   = fractionFlow->getSStar();
  fpwStar = fractionFlow->getFpwStar();
  dtau    = fpwStar / (nPts_ - 2);
  bl_1d    = new double[nPts_];
  tau      = new double[nPts_];
  bl_slope = new double[nPts_ - 1]; 
  for(int i = 0; i < nPts_; ++i) tau[i] = double (i * dtau);
}
double BLSolution:: getFracFlow1(double tau){
   double ss = 0;
   eps = 1.0E-12;
        if(tau <= 0       ) ss = 1.0 - Sor_;          // inlet
   else if(tau == fpwStar ) ss = sStar;               // front
   else if(tau >  fpwStar ) ss = Swc_;                // ahead front
   else                     ss = root( tau, eps );
   return fractionFlow->fw(ss);
}

void BLSolution::cal_BL_solution() {
  for(int i = 0; i < nPts_; ++i) {
          if(tau[i] <= 0       ) bl_1d[i] = 1.0 - Sor_;          // inlet
     else if(tau[i] == fpwStar ) bl_1d[i] = sStar;               // front
     else if(tau[i] >  fpwStar ) bl_1d[i] = Swc_;                // ahead front
     else                        bl_1d[i] = root( tau[i], eps ); // behind front
  }
  for(int i = 0; i < nPts_ - 1; ++i) {
     bl_slope[i] = (bl_1d[i + 1] - bl_1d[i]) / dtau;
  }
}

double BLSolution::getBLSolution(double tau_t) {
   if(tau_t > fpwStar  ) return Swc_;
   if(tau_t == fpwStar ) return sStar;
   int i = (int) (tau_t / dtau);
   if( tau_t >= tau[i] && tau_t <=  tau[i+1] ) {
       return (tau_t - tau[i]) * bl_slope[i] + bl_1d[i];   
   } else {
     cout << tau[i-1] <<' '<<tau[i]  <<' ';
     cout << tau_t    <<' '<<tau[i+1]<<' '<<tau[i+2]<<endl;
     exit(0);
   }
}

// --- calculate the root of (fpw - a = 0) ---
// --- a should between 0 and fpwStar, else no right solution ---
double BLSolution::root(double a, double eps) {
  double x0 = 1. - Sor_, x1 = sStar, xm;
  while( true ) {
    xm = (x0 + x1)/2.0;
    (fractionFlow->fpw(xm) - a > 0)? x1 = xm : x0 = xm;
    if(fabs(x1 - x0) < eps) break;
  }
  return xm;
}

// --- output the BLSolution moments to file for visualization -----
void BLSolution::print(){
  ofstream os1("BL_Solution.out", ios::out);
  for(int i = 0; i < nPts_; ++i) os1 << tau[i]<<' ' << bl_1d[i] << endl;
  os1.close();
}

/*
 // --- solver for BLSolution at time t---------
void BLSolution::get_BL_Table() {
  double startT, endT;
  startT = clock();
  cal_BL_solution(nPts_, tau, bl_1d);
  endT   = clock();
  cerr << "Building BL Table time: " 
       <<(endT-startT)/CLOCKS_PER_SEC << " sec" << endl;
}

void BLSolution::get_BL_Table(int nPts_out, double* tau_out, double* sat_out) {
  double startT, endT;
  startT = clock();
  cal_BL_solution(nPts_out, tau_out, sat_out);
  endT   = clock();
  cerr << "Building BL Table time: " 
       <<(endT-startT)/CLOCKS_PER_SEC << " sec" << endl;
}
*/

// --- calculate deterministic solution ---
// --- tau may be from external ---
/*
void BLSolution::cal_BL_solution(int nPts, double* tau, double* sat) {
  for(int i = 0; i < nPts; ++i) {
          if(tau[i] <= 0            ) sat[i] = 1.0 - Sor_;             // inlet
     else if(tau[i] == fpwStar * t_ ) sat[i] = sStar;                  // front
     else if(tau[i] >  fpwStar * t_ ) sat[i] = Swc_;                   // ahead front
     else                             sat[i] = root( tau[i]/t_, eps ); // behind front
  }
}
*/

