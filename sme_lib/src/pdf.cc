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
#include <math.h>
#include <iostream>

#ifndef _PDFGAUSS_H_
#define _PDFGAUSS_H_

//#define two_pi 3.1415926535897932 * 2.0
double pdfGauss(double& x, double& avg, double& var) {
  double two_pi = 3.1415926535897932 * 2.0;
  return 1./sqrt(two_pi * var) * exp( -0.5 * (x - avg)*(x - avg)/ var);
}

double pdfGauss2(double& x1, double& avg1, double& var1,
		 double& x2, double& avg2, double& var2, double& rho
		) {
  double two_pi = 3.1415926535897932 * 2.0;
  double rho2   = rho * rho;
  double var12  = var1 * var2;
  double tm  = 1./two_pi/sqrt( var12 * (1. - rho2) );
  
  double tm1 = (x1-avg1) * (x1-avg1)/var1 + (x2-avg2) * (x2-avg2)/var2
	     - 2. * rho * (x1-avg1) * (x2-avg2)/sqrt(var12);
  
  return tm * exp( -0.5 /(1. - rho2) * tm1) ;
}

double pdfGauss2(double& x1, double& x2, double& rho) {
  double two_pi = 3.1415926535897932 * 2.0;
  double rho2   = rho * rho;
  double tm  = 1./two_pi/sqrt(1. - rho2);
  double tm1 = (x1 * x1 + x2 * x2 - 2. * rho * x1 * x2);
  return tm * exp( -0.5 /(1. - rho2) * tm1) ;
}


double pdfLognormal(double& x, double& avg, double& var) {
   double two_pi = 3.1415926535897932 * 2.0;
   return 1./sqrt(two_pi * var)/x * exp( -0.5 * (log(x) - avg)*(log(x) - avg)/ var);
}

void transform1(int& len, double* tau, double *var, double* lntau, double* lnvar){
   for(int i = 0; i < len; ++i) {
      lntau[i] = 2. * log(tau[i]) - 0.5 * log(var[i]+ tau[i]*tau[i]);	   
      lnvar[i] = log( var[i]/tau[i]/tau[i] + 1.);
   }
}
void transform2(int& len, double* tau, double *var, double* tau1tau2, double* rho){
		
   int i1i2;
   double std; 
   for(int i2 = 0; i2 < len; ++i2) {
      for(int i1 = 0; i1 < len; ++i1) {
         i1i2 = i1 + i2 * len;
	 if(i1 == i2) {
	    rho[i1i2] = 1.0;
	 }else{
	    std  = sqrt(var[i1] * var[i2]);
	    rho[i1i2] =  log(tau1tau2[i1i2]/tau[i1]/tau[i2] + 1.)/std;
	    /*
	    if(rho[i1i2]!=rho[i1i2]){
		std::cout << "Rho Prob: " <<tau1tau2[i1i2]<<" "<<tau[i1]<<" "<<tau[i2] << " " << std << std::endl;
		exit(1);
	    }
            */
	 }
      }
   }
}

//#undef two_pi
#endif

