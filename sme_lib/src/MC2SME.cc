#include "MC2SME.h"

using namespace std;

MC_TO_SME::MC_TO_SME(int len)
 :length(len){
  initialize();
  read();
  transform1(length, avg, var, avg_log, var_log); 
  transform2(length, avg, var_log, cov, cov_log);
}

MC_TO_SME::~MC_TO_SME(){
  delete [] avg;
  delete [] var;
  delete [] cov;
  delete [] avg_log;
  delete [] var_log;
  delete [] cov_log;
}

void MC_TO_SME::initialize(){
  avg     = new double [length];
  var     = new double [length];
  cov     = new double [length * length];
  avg_log = new double [length];
  var_log = new double [length];
  cov_log = new double [length * length];

  for(int i = 0; i < length; ++i) avg[i] = var[i] = 0.;
  for(int i = 0; i < length * length; ++i) cov[i] = 0.;
  for(int i = 0; i < length; ++i) avg_log[i] = var_log[i] = 0.;
  for(int i = 0; i < length * length; ++i) cov_log[i] = 0.;  
}

void MC_TO_SME::read() {
  char* file="/.automount/matrix/d/a/data/liyl/stochastic/monte_carlo/quarter2d/multi_region/test/Tral_Vari_cell.out";
  ifstream in1;
  in1.open(file);
  if( in1.bad() ) {
      cerr << " Can not open " <<  file  << endl;
      exit(8);
  }

  for(int i = 0; i < length; ++i) {
     in1 >> avg[i] >> var[i];
  }
  in1.close();

  char* file2="/.automount/matrix/d/a/data/liyl/stochastic/monte_carlo/quarter2d/multi_region/test/MCTau_Cov.out";
  ifstream in2;
  in2.open(file2);
  if( in2.bad() ) {
      cerr << " Can not open " <<  file2  << endl;
      exit(8);
  }
  for(int j = 0; j < length; ++j) {
      for(int i = 0; i < length; ++i) {
          in2 >> cov[i + j * length ];
	  //cout<<cov[i + j * length ]<<endl;
      }
  }
  in2.close();
}

void MC_TO_SME::transform1(int& len, double* tau, double *var, double* lntau, double* lnvar){
   for(int i = 0; i < len; ++i) {
      lntau[i] = 2. * log(tau[i]) - 0.5 * log(var[i]+ tau[i]*tau[i]);	   
      lnvar[i] = log( var[i]/tau[i]/tau[i] + 1.);
   }
}

void MC_TO_SME::transform2(int& len, double* tau, double *var, double* tau1tau2, double* rho){
		
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
	 }
      }
   }
}

void MC_TO_SME::print() {
  for(int i = 0; i < length; ++i) {
     cout << avg[i] <<' '<< var[i] <<endl;
  }
}
