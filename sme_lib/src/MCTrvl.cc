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

#include "MCTrvl.h"
#include <stdlib.h>

using namespace std;

MCTrvl::MCTrvl() {
//cout<< "MCTrvl()" << endl;
  NLn = 500;
  Nr  =3500;
  initialize();
  readdata();
  
//test_Trvl_Distribution();
//exit(0);

  adjust_flux();
  calcVeloMoments();
  calcLogTrvlMoments();
  
 // calcTrvlMoments();   
}

MCTrvl::~MCTrvl() {
  delete[] theta;
  delete[] thetaw;
  delete[] tr_rand;
  delete[] v1_rand;
  delete[] v2_rand; 
  delete[] tr_avg;
  delete[] tr_var;
  delete[] qt_avg;
  delete[] r_tau1_tau2;
  delete[] qt2;
}

void MCTrvl::initialize() {
   theta  = new double[NLn];
   thetaw = new double[NLn];
   tr_rand= new double[NLn * Nr];
   v1_rand= new double[NLn * Nr];
   v2_rand= new double[NLn * Nr];
   tr_avg = new double[NLn];
   tr_var = new double[NLn];
   qt_avg = new double[NLn];
   r_tau1_tau2 = new double[NLn * NLn];
   qt2         = new double[NLn * NLn]; 
}

void MCTrvl::readdata() {
   ifstream is_flow;
   is_flow.open("prod_time.tmp",ios::in);
   for(int k = 0; k < Nr; ++k) {
      for(int i = 0; i < NLn; ++i) {
         is_flow>>tr_rand[i + NLn * k]>>v1_rand[i + NLn * k]>>v2_rand[i + NLn * k];
      }
   }
   for(int i = 0; i < NLn; ++i) {
      is_flow>>theta[i]>>thetaw[i];
   }
   is_flow>>r_given;
   is_flow.close();
   cout << "r_given = " << r_given<<endl;
}

void MCTrvl::adjust_flux() {
   double total_flux, tm, factor;
   for(int k = 0; k < Nr; ++k) {
      total_flux = 0.;
      for(int i = 0; i < NLn; ++i) {
         tm = v1_rand[i + k * NLn] * sin( theta[i] ) 
            + v2_rand[i + k * NLn] * cos( theta[i] ); 
          total_flux += tm * thetaw[i] * r_given;  
      }
      factor = 1./total_flux;
      for(int i = 0; i < NLn; ++i) {
         v1_rand[i + k * NLn] *= factor;
         v2_rand[i + k * NLn] *= factor;         
      }
   }
}

void MCTrvl::calcVeloMoments(){
   double* rand_v1= new double[Nr ]; 
   double* rand_v2= new double[Nr ];
   double* q1_avg = new double[NLn];
   double* q1_var = new double[NLn];
   double* q2_avg = new double[NLn];
   double* q2_var = new double[NLn];
   double* th1    = new double[NLn];
   double* th2    = new double[NLn];

   for(int i = 0; i < NLn; ++i) {
      for(int k = 0; k < Nr; ++k) {
          rand_v1[k] = v1_rand[i + k * NLn];
      }
      calcStatistics(Nr, rand_v1, q1_avg[i], q1_var[i]);
      for(int k = 0; k < Nr; ++k) {
          rand_v2[k] = v2_rand[i + k * NLn];
      }
      calcStatistics(Nr, rand_v2, q2_avg[i], q2_var[i]);
      th1[i] = sin( theta[i] );
      th2[i] = cos( theta[i] );
      qt_avg[i]=(q1_avg[i]*th1[i] 
	       + q2_avg[i]*th2[i])*thetaw[i]*r_given;
   }

   double sum1, sum2, sum3, r2 = r_given * r_given;
   double *tmp1 = new double[NLn * NLn];
   double *tmp2 = new double[NLn * NLn];
   double *tmp3 = new double[NLn * NLn]; 
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int k = 0; k < Nr; ++k) rand_v1[k] = v1_rand[j1 + k * NLn];
      for(int j2 = 0; j2 < NLn; ++j2) {
        sum3 = 0.;
         for(int k = 0; k < Nr; ++k) {
            sum3 += rand_v1[k] * v2_rand[j2 + k * NLn];
         }
        tmp3[j2 + j1 * NLn] = sum3/Nr;
      }
   }
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int k = 0; k < Nr; ++k) {
         rand_v1[k] = v1_rand[j1 + k * NLn];
        rand_v2[k] = v2_rand[j1 + k * NLn];
      }
      for(int j2 = 0; j2 <= j1; ++j2) {
        sum1 = 0., sum2 = 0.;
         for(int k = 0; k < Nr; ++k) {
            sum1 += rand_v1[k] * v1_rand[j2 + k * NLn];
            sum2 += rand_v2[k] * v2_rand[j2 + k * NLn];
         }
        tmp1[j2 + j1 * NLn] = sum1/Nr;
        tmp2[j2 + j1 * NLn] = sum2/Nr;
      }
   }
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int j2 = j1 + 1; j2 < NLn; ++j2) {
        tmp1[j2 + j1 * NLn] = tmp1[j1 + j2 * NLn];
        tmp2[j2 + j1 * NLn] = tmp2[j1 + j2 * NLn];
      }
   }

   // calculater qt2
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int j2 = 0; j2 < NLn; ++j2) {
          qt2[j2 + j1 * NLn] =(tmp1[j2 + j1 * NLn] * th1[j1] * th1[j2]
                             + tmp2[j2 + j1 * NLn] * th2[j1] * th2[j2]
                             + tmp3[j2 + j1 * NLn] * th1[j1] * th2[j2] 
                             + tmp3[j1 + j2 * NLn] * th1[j2] * th2[j1]
                               ) * thetaw[j1] * thetaw[j2] * r2;
      }
   }
   delete[] tmp1;
   delete[] tmp2;
   delete[] tmp3;
   delete[] th1;
   delete[] th2;
   delete[] q1_avg;
   delete[] q1_var;
   delete[] q2_avg;
   delete[] q2_var;
   delete[] rand_v1;
   delete[] rand_v2;
   //cout << "delete rand_v2" <<endl;
}

void MCTrvl::calcLogTrvlMoments(){
   double* r = new double[Nr];
   for(int i = 0; i < NLn; ++i) {    
      for(int k = 0; k < Nr; ++k) {
          r[k] = log(tr_rand[i + k * NLn]);
      }
      calcStatistics(Nr, r, tr_avg[i], tr_var[i]);
   }
   double sum;
   //aternative (1) 
   /*
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int k = 0; k < Nr; ++k) {
         r[k] = log(tr_rand[j1 + k * NLn]);
      }
      for(int j2 = 0; j2 < NLn; ++j2) {
        sum = 0.;
         for(int k = 0; k < Nr; ++k) {
             sum += (r[k]-tr_avg[j1])*(log(tr_rand[j2 + k * NLn])-tr_avg[j2]);
         }
        r_tau1_tau2[j2+j1*NLn] = sum/sqrt( tr_var[j1] * tr_var[j2] )/Nr;
        cout << r_tau1_tau2[j2+j1*NLn]<<endl;
      }
   }*/
   //alternative (2) == symmetric and faster
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int k = 0; k < Nr; ++k) {
         r[k] = log(tr_rand[j1 + k * NLn]);
      }
      for(int j2 = 0; j2 <= j1; ++j2) {
        sum = 0.;
         for(int k = 0; k < Nr; ++k) {
             sum += (r[k]-tr_avg[j1])*(log(tr_rand[j2 + k * NLn])-tr_avg[j2]);
         }
        r_tau1_tau2[j2+j1*NLn] = sum/sqrt( tr_var[j1] * tr_var[j2] )/Nr;
      }
   }
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int j2 = j1 + 1; j2 < NLn; ++j2) {
         r_tau1_tau2[j2+j1*NLn] = r_tau1_tau2[j1+j2*NLn];
      }
   }
   delete[] r;
   //cout << "delete r" <<endl;
}

void MCTrvl::calcTrvlMoments(){
   double* r = new double[Nr];
   for(int i = 0; i < NLn; ++i) {    
      for(int k = 0; k < Nr; ++k) {
          r[k] = tr_rand[i + k * NLn];
      }
      calcStatistics(Nr, r, tr_avg[i], tr_var[i]);
   }
   return;
   double sum;
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int k = 0; k < Nr; ++k) {
         r[k] = tr_rand[j1 + k * NLn];
      }
      for(int j2 = 0; j2 <= j1; ++j2) {
         sum = 0.;
         for(int k = 0; k < Nr; ++k) {
             sum += (r[k]-tr_avg[j1])*( tr_rand[j2 + k * NLn]-tr_avg[j2]);
         }
         r_tau1_tau2[j2+j1*NLn] = sum/Nr;
      }
   }
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int j2 = j1 + 1; j2 < NLn; ++j2) {
         r_tau1_tau2[j2+j1*NLn] = r_tau1_tau2[j1+j2*NLn];
      }
   }
   delete[] r;
   
   // transform to lognormal moments
   double rho_tmp, avg_tmp, var_tmp;
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int j2 = 0; j2 < NLn; ++j2) {
          normal2lognormal(tr_avg[j1], tr_var[j1],tr_avg[j2], tr_var[j2],
                           r_tau1_tau2[j2+j1*NLn], avg_tmp, var_tmp, rho_tmp);
	  r_tau1_tau2[j2+j1*NLn] = rho_tmp;
      }
   }
   for(int i = 0; i < NLn; ++i) {
       normal2lognormal(tr_avg[i], tr_var[i], avg_tmp, var_tmp);
       tr_avg[i] = avg_tmp;
       tr_var[i] = var_tmp;
   }
}

void MCTrvl::calcStatistics(int& Nr, double * R, double& avg, double & var){
   double sum = 0.;
   for(int i = 0; i < Nr; ++i) {
      sum += R[i];
   }
   avg = sum/Nr;
   sum = 0.;
   for(int i = 0; i < Nr; ++i) {
       sum += (R[i]-avg)*(R[i]-avg);
   }
   var = sum/Nr;
}

// the relationships between <tau> and <ln tau>
//                   between var. of tau and var. of (ln tau)
void MCTrvl::normal2lognormal(double& tauAvg1, double& tauVar1,
                              double& tauAvg2, double& tauVar2, double& ctau1tau2,
            double& ltauAvg1, double& ltauVar1, double& rho12) { 
   double ltauVar2, tmp;
   ltauVar1 = log(1. + tauVar1/tauAvg1/tauAvg1 );   // vari of ln(tau1)
   ltauVar2 = log(1. + tauVar2/tauAvg2/tauAvg2 );
   ltauAvg1 = log( tauAvg1 ) - 0.5 * ltauVar1;      // mean of ln(tau1)
   tmp      = sqrt(ltauVar1  * ltauVar2);
   rho12   = log(1. + ctau1tau2/tauAvg1/tauAvg2)/tmp;
}

void MCTrvl::normal2lognormal(double&  tauAvg, double&  tauVar,
                              double& ltauAvg, double& ltauVar) { 
   ltauVar = log(1. + tauVar/tauAvg/tauAvg );   // vari of ln(tau)
   ltauAvg = log( tauAvg ) - 0.5 * ltauVar;      // mean of ln(tau)
}

//After testing if the distribution of travel time is lognormal.
//we concluded that they are around 5 -20% error.
void MCTrvl::test_Trvl_Distribution(){
   double* r1 = new double[Nr];
   double* r2 = new double[Nr];
   double* tauAvg  = new double[NLn];
   double* tauVar  = new double[NLn];
   double* ltauAvg = new double[NLn];
   double* ltauVar = new double[NLn];
   double* ctau1tau2 = new double[NLn*NLn];
   double  ltauAvg1, ltauVar1,rho12;
   for(int i = 0; i < NLn; ++i) {    
      for(int k = 0; k < Nr; ++k) {
          r1[k] =     tr_rand[i + k * NLn];
          r2[k] = log(tr_rand[i + k * NLn]); 
      }
      calcStatistics(Nr, r1,  tauAvg[i],  tauVar[i]);
      calcStatistics(Nr, r2, ltauAvg[i], ltauVar[i]);
      normal2lognormal(tauAvg[i], tauVar[i], ltauAvg1, ltauVar1);
      cout <<ltauAvg1 <<' '<<ltauAvg[i]<<' ';
      cout <<ltauVar1 <<' '<<ltauVar[i]<<endl;
   }
   cout<< "test rho"<<endl;
   double sum1, sum2;
   for(int j1 = 0; j1 < NLn; ++j1) {
      for(int k = 0; k < Nr; ++k) {
         r1[k] =      tr_rand[j1 + k * NLn];
         r2[k] = log( tr_rand[j1 + k * NLn]);
      }
      for(int j2 = 0; j2 < NLn; ++j2) {
         sum1 = 0.;
         sum2 = 0.;
         for(int k = 0; k < Nr; ++k) {
             sum1 += (r1[k]- tauAvg[j1])*(     tr_rand[j2 + k * NLn] - tauAvg[j2]);
             sum2 += (r2[k]-ltauAvg[j1])*( log(tr_rand[j2 + k * NLn])-ltauAvg[j2]);
         }
         ctau1tau2[j2+j1*NLn]   = sum1/Nr;
         r_tau1_tau2[j2+j1*NLn] = sum2/sqrt( ltauVar[j1] * ltauVar[j2] )/Nr;
         normal2lognormal(tauAvg[j1], tauVar[j1],tauAvg[j2], tauVar[j2],
                          ctau1tau2[j2+j1*NLn], 
                          ltauAvg1, ltauVar1, rho12);
         double error = (rho12-r_tau1_tau2[j2+j1*NLn])/r_tau1_tau2[j2+j1*NLn];
         cout << rho12 << ' '<<r_tau1_tau2[j2+j1*NLn]<<' '<<error<<endl;
      }
   }
   delete[] ctau1tau2;
   delete[] tauAvg;
   delete[] tauVar;
   delete[] ltauAvg;
   delete[] ltauVar;
   delete[] r1;
   delete[] r2;   
}

void MCTrvl::printTrvlTime(){
   ofstream fileOu("MCTrvlTime.out", ios::out);
   if (!fileOu){
       cerr << "\n File does not exist !\n";
       exit(EXIT_FAILURE);
   }
   for(int i = 0; i < NLn; ++i){
       fileOu<< tr_avg[i] << ' '
	     << tr_var[i] << endl; ;
   }
   fileOu.close();
}
