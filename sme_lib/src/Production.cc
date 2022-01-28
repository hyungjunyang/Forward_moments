#include "Production.h"

Production::Production(Unit u, BLSolution* blSoln, int steps, double t, bool debug_)
: unit(u), blSolnPtr(blSoln), num_steps(steps), t_maximum(t), debug(debug_)
{
   if(debug) cout<<"Production::Production()" << endl;
   if(unit == FIELD) {
      convertUnit();
   }
   
   initialize();

}

Production::Production(Unit u, Domain* domain, BLSolution* blSoln, int steps, double t, bool debug_)
:unit(u), domainPtr(domain), blSolnPtr(blSoln), num_steps(steps), t_maximum(t), debug(debug_)
{
  if(debug) {
     cout<<"Production::Production()" << endl;
     cout<<"volume = " << domainPtr->getVolume() << endl;
  }
  if(unit == FIELD) {
      convertUnit();
  }

  initialize();
}

Production::Production(Unit u, Domain* domain, Flow* flow, BLSolution* blSoln, int steps, double t, bool debug_)
:unit(u), domainPtr(domain), flowPtr(flow), blSolnPtr(blSoln), num_steps(steps), t_maximum(t), debug(debug_)
{
  if(debug) {
     cout<<"Production::Production()" << endl;
     cout<<"volume = " << domainPtr->getVolume() << endl;
  }
  if(unit == FIELD) {
      convertUnit();
  }

  initialize();
}

Production::~Production() {
  if(debug) cout<<"Production::~Production()" << endl;
  //numerical integration
  delete [] xx;
  delete [] xw;
  delete [] gx;
  delete [] gw;
  delete [] xwxw;

  delete [] tau;
  delete [] tau_weight;
  delete [] fw;
  delete [] tau1;
  delete [] tau_weight1;
  delete [] fw1;

  delete [] time_step;
  delete [] water_production;
  delete [] w_cut_avg;
  delete [] w_cut_var;
  delete [] w_cut_var2;
  delete [] o_cut_avg;

  delete [] w_cum_avg;
  delete [] w_cum_var;
  delete [] o_cum_avg;
  delete [] o_cum_var;
  delete [] w_cut_var_1;
  delete [] w_cut_var2_1;

  //water production
  if(wParticlesPtr!=NULL) delete[] wParticlesPtr;
  if(pointPtr!=NULL) delete[] pointPtr;
  if(bTrackPtr!=NULL) delete bTrackPtr;
  if(covYLnTau!=NULL){
    int nWell = prodWellMap.size();
    for(int ii=0;ii<nWell;ii++){
      if(covYLnTau[ii]!=NULL){
	delete[] covYLnTau[ii];
      }
    }
    delete[] covYLnTau;
  }
}

void Production::initialize() {
   num_intg = 50;
   multiple = 5.;
   xx   = new double[num_intg];
   xw   = new double[num_intg];
   gx   = new double[num_intg];
   gw   = new double[num_intg];
   xwxw = new double[num_intg * num_intg];
   double a = -1.;
   double b =  1.;
   double c =  0.;
   gauleg(a, b, xx, xw, num_intg);
   for(int i = 0; i < num_intg; ++i) {
       gx[i] = xx[i] * multiple; 
       gw[i] = xw[i] * multiple * pdfGauss(gx[i], c, b);
       for(int j = 0; j < num_intg; ++j) {
          xwxw[j + i * num_intg] = xw[i] * xw[j] * multiple * multiple;
       }
   }
   
   tau         = new double [num_intg ];
   tau_weight  = new double [num_intg ];
   fw          = new double [num_intg ];
   tau1        = new double [num_intg ];
   tau_weight1 = new double [num_intg ];
   fw1         = new double [num_intg ];
   
   time_step        = new double [num_steps];
   water_production = new double [num_steps];

   w_cut_avg = new double [num_steps];
   w_cut_var = new double [num_steps];
   w_cut_var2= new double [num_steps * num_steps];
   o_cut_avg = new double [num_steps];

   w_cum_avg = new double [num_steps];
   w_cum_var = new double [num_steps];
   o_cum_avg = new double [num_steps];
   o_cum_var = new double [num_steps];
   w_cut_var_1  = new double [num_steps];
   w_cut_var2_1 = new double [num_steps * num_steps];

   time_discretize();

   //water production (Pipat)
   gridPtr = domainPtr->getGrid();
   nWaterProd = nParticles = 0;
   wParticlesPtr = NULL;
   pointPtr = NULL;
   bTrackPtr = NULL;
   covYLnTau = NULL;
}

void Production::time_discretize() {
   dt = t_maximum / (num_steps);
   for(int i = 0; i < num_steps; ++i) {
       time_step[i] = t_maximum * (i+1) / (num_steps);
   }
}

void Production::calcWCumMoments(int NLn, double* tr_avg, double* tr_var,
                                 double* qt_avg, double* cqq, double* rho) { 

   //note: check the algorithm of rho, some symmetric issues 
   debug = false;
   if(debug) {
      ofstream os1("travel_time.out", ios::out);
      for(int i = 0; i < num_steps; ++i) {
         os1 << i <<' '
             << tr_avg[i]<<' ' 
             << tr_var[i]<<' '
             << qt_avg[i]<<endl;
      }
      os1.close();
      ofstream os2("travel_corr.out", ios::out);           
      for(int j2 = 0; j2 < NLn; ++j2) {
          for(int j1 = 0; j1 < NLn; ++j1) {
              os2<< j1 << ' ' << j2 <<' '
                 << rho[j1 + j2 * NLn] << ' '
                 << rho[j2 + j1 * NLn] << ' '
                 << rho[j1 + j2 * NLn] - rho[j2 + j1 * NLn]
                 << endl;
          }
          os2 << endl;
      }
      exit(0);
   }

   double rho_1   = 0.9999;
   double fpwStar = blSolnPtr->getFpwStar();

   double * fFlow    = new double[num_steps * NLn * num_intg];
   double * fFlow_avg= new double[num_steps * NLn];
   double * tau_limit= new double[num_steps];
   double * tt_end   = new double[num_steps * NLn];
   double * tt_bgn   = new double[NLn ];
   double * pdfArray = new double[NLn * num_intg * NLn * num_intg];
   int * blRange2    = new int[num_steps * NLn];
   
   for(int i = 0; i < num_steps; ++i){
       tau_limit[i] = fpwStar * time_step[i]; 
   }
   
   for(int l = 0; l < NLn; l++) {
       tt_bgn[l] = tr_avg[l] - multiple * sqrt( tr_var[l] );
   }
     
   double tmp, std, satu;
   for(int tstep = 0; tstep < num_steps; ++tstep) {
       for(int l = 0; l < NLn; l++) {
           std = sqrt( tr_var[l] );
           tt_end[l + tstep * NLn] = tr_avg[l] + multiple * std;
           if(exp(tt_end[l + tstep * NLn]) > tau_limit[tstep] ) {
                  tt_end[l + tstep * NLn] = log(tau_limit[tstep]);
           }
           blRange2[l + tstep * NLn] = 0; 
           if(tt_bgn[l] < tt_end[l + tstep * NLn] ) {
              for(int kt = 0; kt < num_intg; ++kt) {
                  tmp  = tr_avg[l] + std * gx[kt];
                  satu = blSolnPtr->getBLSolution( exp(tmp)/time_step[tstep] );
                  fFlow[kt + num_intg * l + num_intg * NLn * tstep] = blSolnPtr->getFracFlow(satu);
                  if(fFlow[kt + num_intg * l + num_intg * NLn * tstep] >= 1.e-10) {
                     blRange2[l + tstep * NLn] = kt + 1; 
                  }  
              }
           } 
       }
   }
   
   // mean values
   double sum, fracFlow;
   double total_flow = 0.;
   for(int l = 0; l < NLn; ++l) {
       total_flow += qt_avg[l];
   }
   if(unit == FIELD ) {
      cout << endl;
      cout << "total_flow = " << total_flow / stbPday_m3Psec << " stb/day " << endl;
      cout << endl;
   }
   else {
      cout << "Total Production = " <<  total_flow  << endl;
   }
   
   for(int tstep = 0; tstep < num_steps; ++tstep) {
      sum  = 0.;
      for(int l = 0; l < NLn; ++l) {
         fracFlow = 0.;
            for(int kt = 0; kt < blRange2[l + tstep * NLn]; ++kt ) {
                fracFlow += gw[kt] * fFlow[kt +  num_intg * l + num_intg * NLn * tstep];
            }
            sum += qt_avg[l] * fracFlow;
            fFlow_avg[l + NLn * tstep] = fracFlow;
      }
      w_cut_avg[tstep]    = sum;
      o_cut_avg[tstep]    = total_flow - w_cut_avg[tstep];
      if(tstep == 0) {
         w_cum_avg[tstep] = ( w_cut_avg[tstep] + w_cut_avg[tstep] ) / 2.0 * (time_step[tstep] - 0.);
         o_cum_avg[tstep] = ( o_cut_avg[tstep] + o_cut_avg[tstep] ) / 2.0 * (time_step[tstep] - 0.);
      } else {
         w_cum_avg[tstep] =   w_cum_avg[tstep - 1] 
                          + ( w_cut_avg[tstep - 1] + w_cut_avg[tstep] ) / 2.0
                          * (time_step[tstep] - time_step[tstep - 1]);
         o_cum_avg[tstep] =   o_cum_avg[tstep - 1] 
                          + ( o_cut_avg[tstep - 1] + o_cut_avg[tstep] ) / 2.0 
                          * (time_step[tstep] - time_step[tstep - 1]);
      }
   }

   // pdf calculations
   /*
   for(int j2 = 0; j2 < NLn; ++j2) {
       for(int j1 = 0; j1 < NLn; ++j1) {
           cout << j1 <<" " 
                << j2 <<" " <<  rho[j2 + j1 * NLn] << endl;
       }
   }
   exit(0);
   */

   for(int j2 = 0; j2 < NLn; ++j2) {
       for(int j1 = 0; j1 < NLn; ++j1) {
           if(fabs(rho[j2 + j1 * NLn]) > rho_1) {
              if(j1 != j2) {
                 cout <<"Warning:" << j1 <<" != " <<j2 << ' ' <<"rho = " << rho[j2 + j1 * NLn] <<endl;
              }
              for(int kt = 0; kt < num_intg; ++kt) {
                  for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                      if(kt1 == kt) {
                         pdfArray[kt1 + num_intg * j2 + (kt + num_intg * j1) * num_intg * NLn]
                            = gw[kt] - gw[kt] * gw[kt];
                      } else {
                         pdfArray[kt1 + num_intg * j2 + (kt + num_intg * j1) * num_intg * NLn]
                            = 0 - gw[kt] * gw[kt1];
                      }
                  }
              }
           } else {
              for(int kt  = 0; kt  < num_intg; ++kt ) {
                  for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                      pdfArray[kt1 + num_intg * j2 + (kt + num_intg * j1) * num_intg * NLn]
                      = xwxw[kt1 + kt * num_intg] * pdfGauss2(gx[kt], gx[kt1], rho[j2 + j1 * NLn])
                      - gw[kt] * gw[kt1] ; 
                  }
              }
           }
       }
   }

   // Variance and Covariances of Water and Oil Cuts
   int local_algorithm = 3;
   double sum1;
   if(local_algorithm == 1) {
      // algorithm 1)
      if(debug) cout << "Algorith 1" << endl;
      double * pdfArray1= new double[NLn * num_intg * NLn * num_intg];
      for(int j2 = 0; j2 < NLn; ++j2) {
          for(int j1 = 0; j1 < NLn; ++j1) {
              if(fabs(rho[j2 + j1 * NLn]) > rho_1) {
                 for(int kt = 0; kt < num_intg; ++kt) {
                     for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                         if(kt1 == kt) {
                            pdfArray1[kt1 + num_intg * j2 + (kt + num_intg * j1) * num_intg * NLn]
                               = gw[kt];  
                         } else {
                            pdfArray1[kt1 + num_intg * j2 + (kt + num_intg * j1) * num_intg * NLn]
                               = 0 ; 
                         }
                     }
                 }
              } else {
                 for(int kt  = 0; kt  < num_intg; ++kt ) {
                     for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                         pdfArray1[kt1 + num_intg * j2 + (kt + num_intg * j1) * num_intg * NLn]
                         = xwxw[kt1 + kt * num_intg] * pdfGauss2(gx[kt], gx[kt1], rho[j2 + j1 * NLn]);
                     }
                 }
              }
          }
      }
      
      for(int tstep2 = 0; tstep2 < num_steps; ++tstep2) {
          //cout << "Algo_1 tstep = " << tstep2 << endl;
          for(int tstep1 = 0; tstep1 < num_steps; ++tstep1) {
              sum = 0.;
              sum1= 0.;
              for(int j2 = 0; j2 < NLn; ++j2) {
                  for(int j1 = 0; j1 < NLn; ++j1) {
                      fracFlow = 0.;
                      double fracFlow1= 0.;
                      for(int kt  = 0; kt  < blRange2[j1 + tstep1 * NLn]; ++kt ) {
                      for(int kt1 = 0; kt1 < blRange2[j2 + tstep2 * NLn]; ++kt1) {
                          fracFlow += fFlow[kt +  num_intg * j1 + num_intg * NLn * tstep1]
                                    * fFlow[kt1+  num_intg * j2 + num_intg * NLn * tstep2]
                                    * pdfArray[ kt1 + num_intg * j2 
                                            +  (kt  + num_intg * j1) * num_intg * NLn];
                          fracFlow1+= fFlow[kt +  num_intg * j1 + num_intg * NLn * tstep1]
                                    * fFlow[kt1+  num_intg * j2 + num_intg * NLn * tstep2]
                                    * pdfArray1[ kt1 + num_intg * j2 
                                            +  (kt  + num_intg * j1) * num_intg * NLn];
                      }}
                      //sum  += qt_avg[j1] * qt_avg[j2]  * fracFlow
                      //      + cqq[j1 + j2 * NLn] * fracFlow1;
                      sum  += qt_avg[j1] * qt_avg[j2]  * fracFlow;
                      sum1 += qt_avg[j1] * qt_avg[j2]  * fracFlow;
                  }
              }
              w_cut_var2[tstep1 + tstep2 * num_steps]   = sum;
              w_cut_var2_1[tstep1 + tstep2 * num_steps] = sum1;
          }
      }
      delete [] pdfArray1;
   }else if(local_algorithm == 2) {
      // algorithm 2)
      if(debug) cout << "Algorith 2" << endl;
      for(int tstep2 = 0; tstep2 < num_steps; ++tstep2) {
          //cout << "Algo_2 tstep = " << tstep2 << endl;
          for(int tstep1 = 0; tstep1 < num_steps; ++tstep1) {
              sum = 0.;
              sum1= 0.;
              for(int j2 = 0; j2 < NLn; ++j2) {
                  for(int j1 = 0; j1 < NLn; ++j1) {
                      fracFlow = 0.;
                      for(int kt  = 0; kt  < blRange2[j1 + tstep1 * NLn]; ++kt ) {
                      for(int kt1 = 0; kt1 < blRange2[j2 + tstep2 * NLn]; ++kt1) {
                          fracFlow += fFlow[kt +  num_intg * j1 + num_intg * NLn * tstep1]
                                    * fFlow[kt1+  num_intg * j2 + num_intg * NLn * tstep2]
                                    * pdfArray[ kt1 + num_intg * j2 
                                            +  (kt  + num_intg * j1) * num_intg * NLn];
                      }}
                      
                      /*sum += ( (qt_avg[j1] * qt_avg[j2] + cqq[j1 + j2 * NLn]) * fracFlow
                              + cqq[j1 + j2 * NLn] * fFlow_avg[j1 + NLn * tstep1]
                                                   * fFlow_avg[j2 + NLn * tstep2]
                             );
                      */
                      sum  +=  qt_avg[j1] * qt_avg[j2] * fracFlow;
                      sum1 +=  qt_avg[j1] * qt_avg[j2] * fracFlow;
                  }
              }
              w_cut_var2[tstep1 + tstep2 * num_steps] = sum;
              w_cut_var2_1[tstep1 + tstep2 * num_steps] = sum1;
          }
      }
   } else if(local_algorithm == 3) {
      // algorithm 3)
      if(debug) cout << "Algorith 3" << endl;
      for(int tstep2 = 0; tstep2 < num_steps; ++tstep2) {
          //cout << "Algo_3 tstep = " << tstep2 << endl;
          for(int tstep1 = 0; tstep1 <= tstep2; ++tstep1) {
              sum = 0.;
              sum1= 0.;
              for(int j2 = 0; j2 < NLn; ++j2) {
                  for(int j1 = 0; j1 < NLn; ++j1) {
                      fracFlow = 0.;
                      for(int kt  = 0; kt  < blRange2[j1 + tstep1 * NLn]; ++kt ) {
                      for(int kt1 = 0; kt1 < blRange2[j2 + tstep2 * NLn]; ++kt1) {
                          fracFlow += fFlow[kt +  num_intg * j1 + num_intg * NLn * tstep1]
                                    * fFlow[kt1+  num_intg * j2 + num_intg * NLn * tstep2]
                                    * pdfArray[ kt1 + num_intg * j2 
                                            +  (kt  + num_intg * j1) * num_intg * NLn];
                      }}
                      /*
                      sum += ( (qt_avg[j1] * qt_avg[j2] + cqq[j1 + j2 * NLn]) * fracFlow
                              + cqq[j1 + j2 * NLn] * fFlow_avg[j1 + NLn * tstep1]
                                                   * fFlow_avg[j2 + NLn * tstep2]
                             );*/
                      sum  +=  qt_avg[j1] * qt_avg[j2] * fracFlow;
                      sum1 +=  qt_avg[j1] * qt_avg[j2] * fracFlow;
                  }
              }
              w_cut_var2[tstep1 + tstep2 * num_steps] = sum;
              w_cut_var2[tstep2 + tstep1 * num_steps] = sum;
              w_cut_var2_1[tstep1 + tstep2 * num_steps] = sum1;
              w_cut_var2_1[tstep2 + tstep1 * num_steps] = sum1;
          }
      }
   } else if(local_algorithm == 4) {
      // algorithm 4)
      if(debug) cout << "Algorith 4" << endl;
      for(int tstep2 = 0; tstep2 < num_steps; ++tstep2) {
          //cout << "Algo_4 tstep = " << tstep2 << endl;     
          for(int tstep1 = 0; tstep1 <= tstep2; ++tstep1) {
              sum = 0.;
              sum1= 0.;
              for(int j2 = 0; j2 < NLn; ++j2) {
                  for(int j1 = 0; j1 < j2; ++j1) {
                      fracFlow = 0.;
                      for(int kt  = 0; kt  < blRange2[j1 + tstep1 * NLn]; ++kt ) {
                      for(int kt1 = 0; kt1 < blRange2[j2 + tstep2 * NLn]; ++kt1) {
                          fracFlow += fFlow[kt +  num_intg * j1 + num_intg * NLn * tstep1]
                                    * fFlow[kt1+  num_intg * j2 + num_intg * NLn * tstep2]
                                    * pdfArray[ kt1 + num_intg * j2 
                                            +  (kt  + num_intg * j1) * num_intg * NLn];
                      }}
                      /*
                      sum += ( (qt_avg[j1] * qt_avg[j2] + cqq[j1 + j2 * NLn]) * fracFlow
                              + cqq[j1 + j2 * NLn] * fFlow_avg[j1 + NLn * tstep1]
                                                   * fFlow_avg[j2 + NLn * tstep2]
                             ) * 2.;
                      */
                      sum  += ( (qt_avg[j1] * qt_avg[j2] ) * fracFlow
                              ) * 2.;
                      sum1 += ( (qt_avg[j1] * qt_avg[j2] ) * fracFlow
                              ) * 2.;
                  }
              }
              for(int j2 = 0; j2 < NLn; ++j2) {
                  fracFlow = 0.;
                  for(int kt  = 0; kt  < blRange2[j2 + tstep1 * NLn]; ++kt ) {
                      for(int kt1 = 0; kt1 < blRange2[j2 + tstep2 * NLn]; ++kt1) {
                          fracFlow += fFlow[kt +  num_intg * j2 + num_intg * NLn * tstep1]
                                    * fFlow[kt1+  num_intg * j2 + num_intg * NLn * tstep2]
                                    * pdfArray[ kt1 + num_intg * j2 
                                            +  (kt  + num_intg * j2) * num_intg * NLn];
                      }
                  }
                  /*
                  sum += ( (qt_avg[j2] * qt_avg[j2] + cqq[j2 + j2 * NLn]) * fracFlow
                           + cqq[j2 + j2 * NLn] * fFlow_avg[j2 + NLn * tstep1]
                                                * fFlow_avg[j2 + NLn * tstep2]
                         );
                  */
                  sum  += ( (qt_avg[j2] * qt_avg[j2] ) * fracFlow
                          );

                  sum1 += ( (qt_avg[j2] * qt_avg[j2] ) * fracFlow
                          );
              }
              w_cut_var2[tstep1 + tstep2 * num_steps] = sum;
              w_cut_var2[tstep2 + tstep1 * num_steps] = sum;
              w_cut_var2_1[tstep1 + tstep2 * num_steps] = sum1;
              w_cut_var2_1[tstep2 + tstep1 * num_steps] = sum1;
          }
      }
   } else {
      cout << "No Implementation!!"<<endl;
      exit(0);
   }
   
   // Variance of Cumulative Flux
   double dt2d4 = dt * dt /4.;
   for(int tstep = 0; tstep < num_steps; ++tstep) {
       w_cut_var[tstep]   = w_cut_var2[tstep + tstep * num_steps]; 
       w_cut_var_1[tstep] = w_cut_var2_1[tstep + tstep * num_steps];
       if(tstep == 0) {
          w_cum_var[tstep]  = 0.;
       } else {
          double tmp1 = ( w_cut_var2[tstep     +  tstep      * num_steps] 
                        + w_cut_var2[tstep - 1 + (tstep -1 ) * num_steps]
                        + w_cut_var2[tstep     + (tstep -1 ) * num_steps]
                        + w_cut_var2[tstep - 1 + (tstep    ) * num_steps]
                        );
          for(int tstep1 = 0; tstep1 < tstep - 2; ++tstep1) {
              tmp1 = tmp1 + (  w_cut_var2[tstep1     +  tstep      * num_steps]
                             + w_cut_var2[tstep1 + 1 +  tstep      * num_steps]
                             + w_cut_var2[tstep1     + (tstep - 1) * num_steps]
                             + w_cut_var2[tstep1 + 1 + (tstep - 1) * num_steps]
                             + w_cut_var2[tstep     +  tstep1      * num_steps]
                             + w_cut_var2[tstep     + (tstep1 + 1) * num_steps]
                             + w_cut_var2[tstep - 1 +  tstep1      * num_steps]
                             + w_cut_var2[tstep - 1 + (tstep1 + 1) * num_steps]
                            );
          }
          w_cum_var[tstep] = w_cum_var[tstep - 1] + tmp1 * dt2d4;
       }
       //cout << tstep << ' ' << w_cum_var[tstep] << endl;
   }
    
   delete [] fFlow;
   delete [] fFlow_avg;
   delete [] tau_limit;
   delete [] tt_end;
   delete [] tt_bgn;
   delete [] pdfArray;
   delete [] blRange2;
   
   ofstream os1("WaterCutMoments.out", ios::out);
   ofstream os2("WaterCutMoments2.out", ios::out);
   ofstream os3("WaterCutMoments3.out", ios::out);
   os1 << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(6);
   os2 << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(6);
   double cumu_inj_flow = total_flow * time_step[0];
   for(int i = 0; i < num_steps; ++i) {
       if(unit == FIELD) {
          double std_w_cut = sqrt(w_cut_var[i] / stbPday_m3Psec / stbPday_m3Psec );
          double std_w_cum = sqrt(w_cum_var[i] /(stbPday_m3Psec * day_sec) /(stbPday_m3Psec * day_sec));
          os3 << time_step[i]   /  day_sec << ' '
              << w_cut_var[i] / stbPday_m3Psec / stbPday_m3Psec
              << endl;
          os2 << time_step[i]   /  day_sec << ' '
              << o_cum_avg[i]   / (stbPday_m3Psec * day_sec) << ' '
              << o_cum_avg[i]   / (stbPday_m3Psec * day_sec) - std_w_cum * 1. << ' '
              << o_cum_avg[i]   / (stbPday_m3Psec * day_sec) - std_w_cum * 2. << ' '
              << o_cum_avg[i]   / (stbPday_m3Psec * day_sec) - std_w_cum * 3. << ' '
              << o_cum_avg[i]   / (stbPday_m3Psec * day_sec) + std_w_cum * 1. << ' '
              << o_cum_avg[i]   / (stbPday_m3Psec * day_sec) + std_w_cum * 2. << ' '
              << o_cum_avg[i]   / (stbPday_m3Psec * day_sec) + std_w_cum * 3. << endl;
          os1 << time_step[i]   /  day_sec << ' '
              << o_cut_avg[i]   /  stbPday_m3Psec << ' '
              << w_cut_avg[i]   /  stbPday_m3Psec << ' '
              << std_w_cut  << ' '
              << o_cum_avg[i]   / (stbPday_m3Psec * day_sec) << ' '
              << w_cum_avg[i]   / (stbPday_m3Psec * day_sec) << ' '
              << std_w_cum   / (stbPday_m3Psec * day_sec) /(stbPday_m3Psec * day_sec);
             // << w_cut_var_1[i] / stbPday_m3Psec / stbPday_m3Psec  << ' '
             // << cumu_inj_flow / domainPtr->getVolume() ;
          os1 << endl;
          if(i > 0) {
                 cumu_inj_flow += total_flow * (time_step[i] - time_step[i-1]);
          }
       } else {
          double std_w_cum = sqrt(w_cum_var[i]) ;               
          os3 << cumu_inj_flow / domainPtr->getVolume() << ' '
              << w_cut_var[i] << endl;
          os2 << cumu_inj_flow / domainPtr->getVolume()<< ' '
              << o_cum_avg[i]    << ' '
              << o_cum_avg[i]   - std_w_cum * 1. << ' '
              << o_cum_avg[i]   - std_w_cum * 2. << ' '
              << o_cum_avg[i]   - std_w_cum * 3. << ' '
              << o_cum_avg[i]   + std_w_cum * 1. << ' '
              << o_cum_avg[i]   + std_w_cum * 2. << ' '
              << o_cum_avg[i]   + std_w_cum * 3. << endl;
               
          os1 << cumu_inj_flow / domainPtr->getVolume()<< ' '
              << o_cut_avg[i]    << ' '
              << w_cut_avg[i]    << ' '
              << w_cut_var[i]    << ' '
              << o_cum_avg[i]    << ' '
              << w_cum_avg[i]    << ' '
              << w_cum_var[i]    << ' '
              << w_cut_var_1[i];
          os1 << endl;
          if(i > 0) {
                 cumu_inj_flow += total_flow * (time_step[i] - time_step[i-1]);
          }
       }
   }
   os1.close();
   os2.close();
   os3.close();
   /*
   for(int tstep2 = 0; tstep2 < num_steps; ++tstep2) {
       cout << "tstep2 = " << tstep2 << endl;
       for(int tstep1 = 0; tstep1 <= tstep2; ++tstep1) {
           sum = 0.;
           for(int j2 = 0; j2 < NLn; ++j2) {
               if(tt_bgn[j2] < tt_end[j2 + tstep2 * NLn] ) {
                  for(int j1 = 0; j1 < NLn; ++j1) {
                      if(tt_bgn[j1] < tt_end[j1 + tstep1 * NLn] ) {
                         if(fabs(rho[j2 + j1 * NLn]) > 0.98) {
                            fracFlow = 0.;
                            for(int kt = 0; kt < num_intg; ++kt) {
                                fracFlow += gw[kt]
                                          * fFlow[kt +  num_intg * j1 + num_intg * NLn * tstep1]
                                          * fFlow[kt +  num_intg * j1 + num_intg * NLn * tstep1];
                            }
                         } else {
                            fracFlow = 0.;
                            for(int kt  = 0; kt  < num_intg; ++kt ) {
                                for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                                    tmp = pdfGauss2(gx[kt], gx[kt1], rho[j2 + j1 * NLn]);
                                    fracFlow += xwxw[kt1 + kt * num_intg]
                                              * fFlow[kt +  num_intg * j1 + num_intg * NLn * tstep1]
                                              * fFlow[kt1+  num_intg * j2 + num_intg * NLn * tstep2]
                                              * tmp;
                                }
                            }
                         }
                         sum += qt_avg[j1] * qt_avg[j2] * fracFlow;
                      }
                  }
               }
           }
           w_cut_var2[tstep1 + tstep2 * num_steps] = sum;
           w_cut_var2[tstep2 + tstep1 * num_steps] = w_cut_var2[tstep1 + tstep2 * num_steps];
       }
   }
   */
}

void Production::calcWaterCutMomentFast(int NLn, double* tr_avg, 
                                        double* tr_var, double* qt_avg, 
                                        double* rho) {
   double fracFlow, satu, std, tt1, tt2, tau_limit;
   double sum, tmp;
   double sum1, satu1;
   double std_multiple; 
   double fpwStar = blSolnPtr->getFpwStar();
   double * fFlow = new double[num_intg * NLn];

   for(int tstep = 0; tstep < num_steps; ++tstep) {
      tau_limit = fpwStar * time_step[tstep]; 
      for(int l = 0; l < num_intg * NLn; l++) fFlow[l] = 0;
          
      sum  = 0.;
      sum1 = 0.;
      for(int j = 0; j < NLn; ++j) {
         tmp   = exp( tr_avg[j] + 0.5 * tr_var[j] );
         satu1 = blSolnPtr->getBLSolution( tmp/time_step[tstep]);
         sum1 += qt_avg[j] * blSolnPtr->getFracFlow(satu1);
              
         std = sqrt( tr_var[j] );
         std_multiple = multiple * std;
         tt1 = tr_avg[j] - std_multiple;
         tt2 = tr_avg[j] + std_multiple;
         if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
         
         fracFlow = 0.;
         if(tt1 < tt2 ) {
            for(int kt = 0; kt < num_intg; ++kt) {
                tmp  = tr_avg[j] + std * gx[kt];
                satu = blSolnPtr->getBLSolution( exp(tmp)/time_step[tstep] );

                fFlow[kt +  num_intg * j] = blSolnPtr->getFracFlow(satu);
                fracFlow += gw[kt] * fFlow[kt + num_intg * j];
            }
            sum += qt_avg[j] * fracFlow;
         } 
      }
      w_cut_avg[tstep]    = sum;
      if(tstep == 0) {
         w_cum_avg[tstep] = 0;
      } else {
         w_cum_avg[tstep] =   w_cum_avg[tstep - 1] 
                          + ( w_cut_avg[tstep - 1] + w_cut_avg[tstep] ) / 2.0 * dt;
      }
      
      water_production[tstep] = sum1;

      //variance
      sum = 0.;
      for(int j1 = 0; j1 < NLn; ++j1) {
         std = sqrt( tr_var[j1] );
         std_multiple = multiple * std;
         tt1 = tr_avg[j1] - std_multiple;
         tt2 = tr_avg[j1] + std_multiple;         
         if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
         
         if(tt1 < tt2 ) {
            for(int j2 = 0; j2 < NLn; ++j2) { 
                std = sqrt( tr_var[j2] );
                std_multiple = multiple * std;
                tt1 = tr_avg[j2] - std_multiple;
                tt2 = tr_avg[j2] + std_multiple;                
                if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
                
                if(tt1 < tt2 ) {
                   if(fabs(rho[j2 + j1 * NLn]) > 0.98) {
                      fracFlow = 0.;
                      for(int kt = 0; kt < num_intg; ++kt) {
                         fracFlow += gw[kt] * fFlow[kt + num_intg * j1]
                                            * fFlow[kt + num_intg * j1];
                      }
                   } else {
                      fracFlow = 0.;
                      for(int kt  = 0; kt  < num_intg; ++kt ) {
                      for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                             tmp = pdfGauss2(gx[kt], gx[kt1], rho[j2 + j1 * NLn]);
                             fracFlow += xwxw[kt1 + kt * num_intg]
                                       * fFlow[kt  + num_intg * j1] 
                                       * fFlow[kt1 + num_intg * j2]
                                       * tmp;
                      }
                      }
                   }
                   sum += qt_avg[j1] * qt_avg[j2] * fracFlow;
                }
            }
         } 
      }
      
      w_cut_var[tstep] = sum - w_cut_avg[tstep] * w_cut_avg[tstep];
      if(debug) cout << tstep <<' ' << w_cut_var[tstep] << endl;
   }

   delete [] fFlow;

   ofstream os1("WaterCutMoments.out", ios::out);
   os1 << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(6);
   for(int i = 0; i < num_steps; ++i) {
       os1 << time_step[i]/domainPtr->getVolume()<< ' '
           << w_cut_avg[i]    << ' '
           << w_cut_var[i]    << ' '
           << w_cum_avg[i]    << ' '
           << w_cum_var[i]    << ' '
           << water_production[i] <<' ';
       os1 << endl;
   }
   os1.close();
}

void Production::calcWaterCutMomentFast(int NLn, double* tr_avg, 
                                    double* tr_var, double* qt_avg) {
   double fracFlow, satu, std, tt1, tt2, tau_limit;
   double sum, tmp;
   double std_multiple; 
   double fpwStar = blSolnPtr->getFpwStar();
   for(int tstep = 0; tstep < num_steps; ++tstep) {
      tau_limit = fpwStar * time_step[tstep]; 
      sum = 0.;
      for(int j = 0; j < NLn; ++j) {
         std = sqrt( tr_var[j] );
         std_multiple = multiple * std;
         tt1 = tr_avg[j] - std_multiple;
         tt2 = tr_avg[j] + std_multiple;
         if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
         fracFlow = 0.;
         if(tt1 < tt2 ) {
            for(int kt = 0; kt < num_intg; ++kt) {
                tmp  = tr_avg[j] + std *  gx[kt];
                satu = blSolnPtr->getBLSolution( exp(tmp)/time_step[tstep] );
                fracFlow += gw[kt] * blSolnPtr->getFracFlow(satu);
            }
            sum += qt_avg[j] * fracFlow;
         } 
      }
      w_cut_avg[tstep] = sum;
   }
 
   ofstream os1("WaterCutMean.out", ios::out);
   for(int i = 0; i < num_steps; ++i) {
      os1 << time_step[i]/domainPtr->getVolume()<<' '
          << w_cut_avg[i]  <<endl;
   }
   os1.close();
}

void Production::calcWaterCutMoments(int NLn, double* tr_avg, 
                 double* tr_var, double* qt_avg) {
   double fracFlow, satu, std, tt1, tt2, tau_limit;
   double sum, tmp;
   double fpwStar = blSolnPtr->getFpwStar();
   for(int tstep = 0; tstep < num_steps; ++tstep) {
      tau_limit = fpwStar * time_step[tstep]; 
      sum = 0.;
      for(int j = 0; j < NLn; ++j) {
         std = sqrt( tr_var[j] );
         tt1 = tr_avg[j] - 5. * std;
         tt2 = tr_avg[j] + 5. * std; 
         if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
         fracFlow = 0.;
         if(tt1 < tt2 ) {
            gauleg(tt1, tt2, tau, tau_weight, num_intg);
            for(int kt = 0; kt < num_intg; ++kt) {
                satu = blSolnPtr->getBLSolution(exp(tau[kt])/time_step[tstep]);
                fracFlow += tau_weight[kt] * blSolnPtr->getFracFlow(satu)
                            * pdfGauss(tau[kt],tr_avg[j],tr_var[j]);
            }
            sum += qt_avg[j] * fracFlow;
         } 
      }
      w_cut_avg[tstep] = sum;
   }
 
   ofstream os1("WaterCutMean.out", ios::out);
   for(int i = 0; i < num_steps; ++i) {
      os1 << time_step[i]/domainPtr->getVolume()<<' '
          << w_cut_avg[i]  <<endl;
   }
   os1.close();
}

void Production::calcWaterCutMoments(int NLn, double* tr_avg, 
                 double* tr_var, double* qt_avg, double* rho) {
   double fracFlow, satu, std, tt1, tt2, tau_limit;
   double sum, tmp;
   double fpwStar = blSolnPtr->getFpwStar();
   for(int tstep = 0; tstep < num_steps; ++tstep) {
      tau_limit = fpwStar * time_step[tstep]; 
      sum = 0.;
      for(int j = 0; j < NLn; ++j) {
         std = sqrt( tr_var[j] );
         tt1 = tr_avg[j] - 5. * std;
         tt2 = tr_avg[j] + 5. * std; 
         if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
         fracFlow = 0.;
         if(tt1 < tt2 ) {
            gauleg(tt1, tt2, tau, tau_weight, num_intg);
            for(int kt = 0; kt < num_intg; ++kt) {
                satu = blSolnPtr->getBLSolution(exp(tau[kt])/time_step[tstep]);
                  fracFlow += tau_weight[kt] * blSolnPtr->getFracFlow(satu)
                            * pdfGauss(tau[kt],tr_avg[j],tr_var[j]);
            }
            sum += qt_avg[j] * fracFlow;
         } 
      }
      w_cut_avg[tstep] = sum;
   }

   // variance
   for(int tstep = 0; tstep < num_steps; ++tstep) {
      tau_limit = fpwStar * time_step[tstep]; 
      sum = 0.;
      for(int j1 = 0; j1 < NLn; ++j1) {
         std = sqrt( tr_var[j1] );
         tt1 = tr_avg[j1] - 5. * std;
         tt2 = tr_avg[j1] + 5. * std; 
         if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
         
         if(tt1 < tt2 ) {
            gauleg(tt1, tt2, tau, tau_weight, num_intg);
            for(int kt = 0; kt < num_intg; ++kt) {
                satu   = blSolnPtr->getBLSolution( exp(tau[kt]) /time_step[tstep] );
                fw[kt] = blSolnPtr->getFracFlow(satu);
            }
            for(int j2 = 0; j2 < NLn; ++j2) { 
                std = sqrt( tr_var[j2] );
                tt1 = tr_avg[j2] - 5. * std;
                tt2 = tr_avg[j2] + 5. * std; 
                if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
                if(tt1 < tt2 ) {
                   gauleg(tt1, tt2, tau1, tau_weight1, num_intg);
                   for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                       satu     = blSolnPtr->getBLSolution(exp(tau1[kt1])/time_step[tstep]);
                       fw1[kt1] = blSolnPtr->getFracFlow(satu);
                   }
                   if(fabs(rho[j2 + j1 * NLn]) > 0.98) {
                      fracFlow = 0.;
                      for(int kt = 0; kt < num_intg; ++kt) {
                         fracFlow += tau_weight[kt] * fw[kt] * fw[kt]
                                   * pdfGauss(tau[kt],tr_avg[j1],tr_var[j1]);
                      }
                   } else {
                      fracFlow = 0.;
                      for(int kt = 0; kt < num_intg; ++kt) {
                         for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                             tmp = pdfGauss2(tau [kt ],tr_avg[j1],tr_var[j1],
                                             tau1[kt1],tr_avg[j2],tr_var[j2],
                                             rho[j2 + j1 * NLn]);
                             fracFlow += tau_weight [kt]  * fw [kt] * 
                                         tau_weight1[kt1] * fw1[kt1]* tmp;
                         }
                      }
                   }
                   sum += qt_avg[j1] * qt_avg[j2] * fracFlow;
                }
            }
         } 
      }
      w_cut_var[tstep] = sum - w_cut_avg[tstep] * w_cut_avg[tstep];
      if(debug) cout << tstep <<' ' << w_cut_var[tstep] << endl;
   }

   ofstream os1("WaterCutMoments.out", ios::out);
   for(int i = 0; i < num_steps; ++i) {
      os1 << time_step[i]/domainPtr->getVolume() << ' '
          << w_cut_avg[i] << ' '
          << w_cut_var[i] <<endl;
   }
   os1.close();
}

void Production::calcWaterCutMoments(int NLn, double* tr_avg, double* tr_var,  
                                     double* qt_avg, double* qt2, double* rho) {
   
   double fracFlow, satu, std, tt1, tt2, tau_limit;
   double sum, tmp;
   double fpwStar = blSolnPtr->getFpwStar();
   for(int tstep = 0; tstep < num_steps; ++tstep) {
      tau_limit = fpwStar * time_step[tstep]; 
      sum = 0.;
      for(int j = 0; j < NLn; ++j) {
         std = sqrt( tr_var[j] );
         tt1 = tr_avg[j] - 5. * std;
         tt2 = tr_avg[j] + 5. * std; 
         if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
         fracFlow = 0.;
         if(tt1 < tt2 ) {
            gauleg(tt1, tt2, tau, tau_weight, num_intg);
            for(int kt = 0; kt < num_intg; ++kt) {
                satu = blSolnPtr->getBLSolution(exp(tau[kt])/time_step[tstep]);
                fracFlow += tau_weight[kt] * blSolnPtr->getFracFlow(satu)
                            * pdfGauss(tau[kt],tr_avg[j],tr_var[j]);
            }
            sum += qt_avg[j] * fracFlow;
         } 
      }
      w_cut_avg[tstep] = sum;
   }
   // variance
   for(int tstep = 0; tstep < num_steps; ++tstep) {
      tau_limit = fpwStar * time_step[tstep]; 
      sum = 0.;
      for(int j1 = 0; j1 < NLn; ++j1) {
         std = sqrt( tr_var[j1] );
         tt1 = tr_avg[j1] - 5. * std;
         tt2 = tr_avg[j1] + 5. * std; 
         if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
         
         if(tt1 < tt2 ) {
            gauleg(tt1, tt2, tau, tau_weight, num_intg);
            for(int kt = 0; kt < num_intg; ++kt) {
                satu   = blSolnPtr->getBLSolution( exp(tau[kt]) /time_step[tstep] );
                fw[kt] = blSolnPtr->getFracFlow(satu);
            }
            for(int j2 = 0; j2 < NLn; ++j2) { 
                std = sqrt( tr_var[j2] );
                tt1 = tr_avg[j2] - 5. * std;
                tt2 = tr_avg[j2] + 5. * std; 
                if( exp(tt2) > tau_limit ) tt2 = log(tau_limit);
                if(tt1 < tt2 ) {
                   gauleg(tt1, tt2, tau1, tau_weight1, num_intg);
                   for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                       satu     = blSolnPtr->getBLSolution(exp(tau1[kt1])/time_step[tstep]);
                       fw1[kt1] = blSolnPtr->getFracFlow(satu);
                   }
                   if(fabs(rho[j2 + j1 * NLn]) > 0.98) {
                      fracFlow = 0.;
                      for(int kt = 0; kt < num_intg; ++kt) {
                         fracFlow += tau_weight[kt] * fw[kt] * fw[kt]
                                   * pdfGauss(tau[kt],tr_avg[j1],tr_var[j1]);
                      }
                   } else {
                      fracFlow = 0.;
                      for(int kt = 0; kt < num_intg; ++kt) {
                         for(int kt1 = 0; kt1 < num_intg; ++kt1) {
                             tmp = pdfGauss2(tau [kt ],tr_avg[j1],tr_var[j1],
                                             tau1[kt1],tr_avg[j2],tr_var[j2],
                                             rho[j2 + j1 * NLn]);
                             fracFlow += tau_weight [kt]  * fw [kt] * 
                                         tau_weight1[kt1] * fw1[kt1]* tmp;
                         }
                      }
                   }
                   sum += qt2[j2 + j1 * NLn] * fracFlow;
                }
            }
         } 
      }
      w_cut_var[tstep] = sum - w_cut_avg[tstep] * w_cut_avg[tstep];
      if(debug) cout << tstep <<' ' << w_cut_var[tstep] << endl;
   }

   ofstream os1("WaterCutMoments.out", ios::out);
   for(int i = 0; i < num_steps; ++i) {
      os1 << time_step[i]/domainPtr->getVolume() << ' '
          << w_cut_avg[i] << ' '
          << w_cut_var[i] <<endl;
   }
   os1.close();
}

void Production::calcWaterCut( ifstream &is, char *Dir ) {
   double satu;
   double *v_total;
   for(int i = 0; i < domainPtr->getNumWells(); ++i) {
       if( domainPtr->getWellArc(i)->isInjcWell() ) {
           int nstrl_well = domainPtr->getWellArc(i)->length();
           pt_start = new Point[nstrl_well];
           for(int j = 0; j < nstrl_well; ++j) {
               pt_start[j] = (*domainPtr->getWellArc(i))[j];
           }
           ptclTrackPtr  = new ParticleTrack(domainPtr, nstrl_well, pt_start, false);
           ptclTrackPtr->goTrack();
           double *tauAvg = new double [nstrl_well];
           for(int j = 0; j < nstrl_well; ++j) {
               tauAvg[j] = ptclTrackPtr->getParticle(j)->getLastTravelTime();
           }
           //v_total = domainPtr->getWellArc(i)->getTotalV();
	   cout <<"v_total is modified, but uncompleted yet!!!" << endl;
	   cout <<"please revisit the code" << endl;
	   exit(0);
           for(int tstep = 0; tstep < num_steps; ++tstep) {
               double sum = 0.;
               for(int j = 0; j < nstrl_well; ++j) {
                   satu = blSolnPtr->getBLSolution(tauAvg[j]/time_step[tstep]);
                   sum += v_total[j] * blSolnPtr->getFracFlow(satu);
               }
               water_production[tstep] = sum;
           }
           delete [] tauAvg;
           delete [] pt_start;
           delete ptclTrackPtr;
       }
   }
   ofstream os1("Water_Cut.out", ios::out);
   for(int i = 0; i < num_steps; ++i) {
      os1 << time_step[i]<<' '<< water_production[i] << endl;
   }
   os1.close();
}

// --- convert field unit to metric unit -------------
void Production::convertUnit() 
{  t_maximum *= day_sec;
}

double Production::CalcTotalProdRateMean(int wInd)
{
	Well* wellsPtr = domainPtr->getWellArr();
	int wLen = wellsPtr[wInd].getWellLen();
	double rateMean = 0;
	for(int ib=0;ib<wLen;ib++)
	{
		int ii = wellsPtr[wInd].getWellIJK(0,ib);
		int jj = wellsPtr[wInd].getWellIJK(1,ib);
		int kk = wellsPtr[wInd].getWellIJK(2,ib);
		double dx = domainPtr->getBlock(ii,jj,kk)->getDx();
		double dy = domainPtr->getBlock(ii,jj,kk)->getDy();
		double dz = domainPtr->getBlock(ii,jj,kk)->getDz();
		double poro = domainPtr->getPoro(ii,jj,kk);
		
		rateMean += domainPtr->getBlock(ii,jj,kk)->getVx1()*poro*dy*dz;
		rateMean -= domainPtr->getBlock(ii,jj,kk)->getVx2()*poro*dy*dz;
		rateMean += domainPtr->getBlock(ii,jj,kk)->getVy1()*poro*dx*dz;
		rateMean -= domainPtr->getBlock(ii,jj,kk)->getVy2()*poro*dx*dz;
	}
	return(rateMean);
}

double Production::CalcTotalProdRateCov(int wInd1, int wInd2)
{
	Well* wellsPtr = domainPtr->getWellArr();
	int wLen1 = wellsPtr[wInd1].getWellLen();
	int wLen2 = wellsPtr[wInd2].getWellLen();
	double* Cv1v1 = domainPtr->getCv1v1();
	double* Cv1v2 = domainPtr->getCv1v2();
	double* Cv2v2 = domainPtr->getCv2v2();
	int nx = domainPtr->getNx();
	int ny = domainPtr->getNy();
	double rateCov =0;
	for(int ib1=0;ib1<wLen1;ib1++)
	{
		int ii1 = wellsPtr[wInd1].getWellIJK(0,ib1);
		int jj1 = wellsPtr[wInd1].getWellIJK(1,ib1);
		int kk1 = wellsPtr[wInd1].getWellIJK(2,ib1);
		int ijk1 = domainPtr->getIndex(ii1,jj1,kk1);
		double dx1 = domainPtr->getBlock(ii1,jj1,kk1)->getDx();
		double dy1 = domainPtr->getBlock(ii1,jj1,kk1)->getDy();
		double dz1 = domainPtr->getBlock(ii1,jj1,kk1)->getDz();
		double poro1 = domainPtr->getPoro(ii1,jj1,kk1);
		for(int ib2=0;ib2<wLen2;ib2++)
		{
			int ii2 = wellsPtr[wInd2].getWellIJK(0,ib2);
			int jj2 = wellsPtr[wInd2].getWellIJK(1,ib2);
			int kk2 = wellsPtr[wInd2].getWellIJK(2,ib2);
			int ijk2 = domainPtr->getIndex(ii2,jj2,kk2);
			double dx2 = domainPtr->getBlock(ii2,jj2,kk2)->getDx();
			double dy2 = domainPtr->getBlock(ii2,jj2,kk2)->getDy();
			double dz2 = domainPtr->getBlock(ii2,jj2,kk2)->getDz();
			double poro2 = domainPtr->getPoro(ii2,jj2,kk2);
			int nNode = domainPtr->getGrid()->getnNode();

			//flux from the top and bottom of the well is not included here
			//This assumes no flux from top and bottom for block 1 to n-2 (need Cv3vi)
			//This assumes no flow boundary condition
			if(ii1>0)
			{
				if(ii2>0){rateCov += Cv1v1[(domainPtr->getIndex(ii1-1,jj1,kk1))+((domainPtr->getIndex(ii2-1,jj2,kk2))*nNode)]*poro1*poro2*dy1*dz1*dy2*dz2;}
				if(ii2<nx-1){rateCov -= Cv1v1[(domainPtr->getIndex(ii1-1,jj1,kk1))+((domainPtr->getIndex(ii2,jj2,kk2))*nNode)]*poro1*poro2*dy1*dz1*dy2*dz2;}
				if(jj2>0){rateCov += Cv1v2[(domainPtr->getIndex(ii1-1,jj1,kk1))+((domainPtr->getIndex(ii2,jj2-1,kk2))*nNode)]*poro1*poro2*dy1*dz1*dx2*dz2;}
				if(jj2<ny-1){rateCov -= Cv1v2[(domainPtr->getIndex(ii1-1,jj1,kk1))+((domainPtr->getIndex(ii2,jj2,kk2))*nNode)]*poro1*poro2*dy1*dz1*dx2*dz2;}
			}
			if(ii1<nx-1)
			{
				if(ii2>0){rateCov -= Cv1v1[(domainPtr->getIndex(ii1,jj1,kk1))+((domainPtr->getIndex(ii2-1,jj2,kk2))*nNode)]*poro1*poro2*dy1*dz1*dy2*dz2;}
				if(ii2<nx-1){rateCov += Cv1v1[(domainPtr->getIndex(ii1,jj1,kk1))+((domainPtr->getIndex(ii2,jj2,kk2))*nNode)]*poro1*poro2*dy1*dz1*dy2*dz2;}
				if(jj2>0){rateCov -= Cv1v2[(domainPtr->getIndex(ii1,jj1,kk1))+((domainPtr->getIndex(ii2,jj2-1,kk2))*nNode)]*poro1*poro2*dy1*dz1*dx2*dz2;}
				if(jj2<ny-1){rateCov += Cv1v2[(domainPtr->getIndex(ii1,jj1,kk1))+((domainPtr->getIndex(ii2,jj2,kk2))*nNode)]*poro1*poro2*dy1*dz1*dx2*dz2;}
			}
			if(jj1>0)
			{
				if(ii2>0){rateCov += Cv1v2[(domainPtr->getIndex(ii2-1,jj2,kk2))+((domainPtr->getIndex(ii1,jj1-1,kk1))*nNode)]*poro1*poro2*dx1*dz1*dy2*dz2;}
				if(ii2<nx-1){rateCov -= Cv1v2[(domainPtr->getIndex(ii2,jj2,kk2))+((domainPtr->getIndex(ii1,jj1-1,kk1))*nNode)]*poro1*poro2*dx1*dz1*dy2*dz2;}
				if(jj2>0){rateCov += Cv2v2[(domainPtr->getIndex(ii1,jj1-1,kk1))+((domainPtr->getIndex(ii2,jj2-1,kk2))*nNode)]*poro1*poro2*dx1*dz1*dx2*dz2;}
				if(jj2<ny-1){rateCov -= Cv2v2[(domainPtr->getIndex(ii1,jj1-1,kk1))+((domainPtr->getIndex(ii2,jj2,kk2))*nNode)]*poro1*poro2*dx1*dz1*dx2*dz2;}
			}
			if(jj1<ny-1)
			{
				if(ii2>0){rateCov -= Cv1v2[(domainPtr->getIndex(ii2-1,jj2,kk2))+((domainPtr->getIndex(ii1,jj1,kk1))*nNode)]*poro1*poro2*dx1*dz1*dy2*dz2;}
				if(ii2<nx-1){rateCov += Cv1v2[(domainPtr->getIndex(ii2,jj2,kk2))+((domainPtr->getIndex(ii1,jj1,kk1))*nNode)]*poro1*poro2*dx1*dz1*dy2*dz2;}
				if(jj2>0){rateCov -= Cv2v2[(domainPtr->getIndex(ii1,jj1,kk1))+((domainPtr->getIndex(ii2,jj2-1,kk2))*nNode)]*poro1*poro2*dx1*dz1*dx2*dz2;}
				if(jj2<ny-1){rateCov += Cv2v2[(domainPtr->getIndex(ii1,jj1,kk1))+((domainPtr->getIndex(ii2,jj2,kk2))*nNode)]*poro1*poro2*dx1*dz1*dx2*dz2;}
			}
		}
	}
	return(rateCov);
}

void Production::ProdWellIs(int wInd){
  map<int,int>::iterator it = prodWellMap.find(wInd);
  if(it==prodWellMap.end()){
    prodWellMap[wInd] = nWaterProd++;
  }
}

void Production::AddAllProdWells(){
  Well* wellsPtr = domainPtr->getWellArr();
  for(int wInd=0;wInd<domainPtr->getNumWells();wInd++){
    if(wellsPtr[wInd].isProd()) ProdWellIs(wInd);
  }
}

void Production::SetParticlesAtProdWells(){
  if(nWaterProd==0) AddAllProdWells();
  int nWell = prodWellMap.size();
  wParticlesPtr = new WellParticles[nWell];
  pointPtr = new Point[nWell*4*int_section_per_side*int_quad_point];
  double* rloc, * weight;
  GaussianQuadrature(int_quad_point,rloc,weight);

  Well* wellsPtr = domainPtr->getWellArr();
  int di[6] = {-1,1,0,0,0,0};
  int dj[6] = {0,0,-1,1,0,0};
  int dk[6] = {0,0,0,0,-1,1};
  for(map<int,int>::iterator itw=prodWellMap.begin();itw!=prodWellMap.end();itw++){
    int wi,wj,wk;
    wi = wellsPtr[itw->first].getWellIJK(0,0);
    wj = wellsPtr[itw->first].getWellIJK(1,0);
    wk = wellsPtr[itw->first].getWellIJK(2,0);
    double wx,wy,wz;
    wx = gridPtr->getX(wi);
    wy = gridPtr->getY(wj);
    wz = gridPtr->getZ(wk);
    double dx = domainPtr->getBlock(wi,wj,wk)->getDx();
    double dy = domainPtr->getBlock(wi,wj,wk)->getDy();
    wParticlesPtr[itw->second].start = nParticles;

    for(int side=0;side<4;side++){
      if(domainPtr->IsInDomain(wi+di[side],wj+dj[side],wk+dk[side])){
	for(int sec=0;sec<int_section_per_side;sec++){
	  for(int qp=0;qp<int_quad_point;qp++){
	    //Calculate position using general relation, works for 2D now
	    Point tempPoint(wx+(di[side]*dx/2)+(abs(dj[side])*(-(dx/2)+(sec*dx/int_section_per_side)+(dx/int_section_per_side/2*(rloc[qp]+1)))),
			    wy+(abs(di[side])*(-(dy/2)+(sec*dy/int_section_per_side)+(dy/int_section_per_side/2*(rloc[qp]+1))))+(dj[side]*dy/2),
			    0,  //To be consistent
			    wi,wj,wk);
	    pointPtr[nParticles++] = tempPoint;
	  }
	}
	wParticlesPtr[itw->second].side[side] = true;
      }
      else wParticlesPtr[itw->second].side[side] = false;
    }
  }

  delete[] rloc;
  delete[] weight;

  //initiate covYLnTau
  covYLnTau = new double*[nWell];
  for(int ii=0;ii<nWell;ii++){
    covYLnTau[ii] = NULL;
  }
}

void Production::LaunchParticles(){
  if(bTrackPtr!=NULL) delete bTrackPtr;

  bTrackPtr = new ParticleTrack(domainPtr,nParticles,pointPtr,true);
  bTrackPtr->goTrackBackward();
  bTrackPtr->calcTrvlTimeMomentsNoWell();
  bTrackPtr->printTravelTime();
}

void Production::CalcCovYLnTau(int wInd){
  SwCondiInverse* swCondPtr = new SwCondiInverse(domainPtr,flowPtr);
  swCondPtr->BTrackPtrIs(bTrackPtr);
  Perm* permPtr = flowPtr->getPermPtr();

  int nStrl = 0;
  for(int ii=0;ii<4;ii++){
    if(wParticlesPtr[prodWellMap[wInd]].side[ii]) nStrl += int_section_per_side*int_quad_point;
  }

  if(covYLnTau[prodWellMap[wInd]]!=NULL){
    delete[] covYLnTau[prodWellMap[wInd]];
  }
  covYLnTau[prodWellMap[wInd]] = new double[nStrl*gridPtr->getNumPermNode()];
  
  int start = wParticlesPtr[prodWellMap[wInd]].start;
  for(int istrl=0;istrl<nStrl;istrl++){
    swCondPtr->CollectDataS(start+istrl);
    for(int iy=0;iy<gridPtr->getNumPermNode();iy++){
      for(int ip=0;ip<bTrackPtr->getParticle(start+istrl)->getLength();ip++){
	swCondPtr->CalcCov_YVxi_YVeta_YEta(start+istrl,iy,ip);
      }
      double covYTau = swCondPtr->CalcCovYTau(start+istrl);
      covYLnTau[prodWellMap[wInd]][istrl*gridPtr->getNumPermNode()+iy] =
	swCondPtr->CalcCovYLnTau(bTrackPtr->getParticle(start+istrl)->getLastTravelTime(),covYTau);
    }
    swCondPtr->DeallocateDataS();
  }

  

  delete swCondPtr;
}

void Production::ClearCovYLnTau(){
  int nWell = prodWellMap.size();
  for(int ii=0;ii<nWell;ii++){
    if(covYLnTau[ii]!=NULL){
      delete[] covYLnTau[ii];
      covYLnTau[ii] = NULL;
    }
  }
}

double Production::CalcWaterProdRateMean(int wInd, double wTime){
  double qwMn = 0;
  Well* wellsPtr = domainPtr->getWellArr();

  //The current implementation allows only one-block well
  int wi = wellsPtr[wInd].getWellIJK(0,0);
  int wj = wellsPtr[wInd].getWellIJK(1,0);
  int wk = wellsPtr[wInd].getWellIJK(2,0);
  double dx = domainPtr->getBlock(wi,wj,wk)->getDx();
  double dy = domainPtr->getBlock(wi,wj,wk)->getDy();
  double dz = domainPtr->getBlock(wi,wj,wk)->getDz();
  double poro = domainPtr->getPoro(wi,wj,wk);
  
  double darcyV[4];
  darcyV[0] = domainPtr->getBlock(wi,wj,wk)->getVx1()*poro;
  darcyV[1] = domainPtr->getBlock(wi,wj,wk)->getVx2()*poro;
  darcyV[2] = domainPtr->getBlock(wi,wj,wk)->getVy1()*poro;
  darcyV[3] = domainPtr->getBlock(wi,wj,wk)->getVy2()*poro;
  double lSection[4] = {dy/int_section_per_side,dy/int_section_per_side,dx/int_section_per_side,dx/int_section_per_side};

  double* rlocw, * weightw;
  GaussianQuadrature(int_quad_point,rlocw,weightw);
  double sDRange = 3.5, intRange = 1;
  int qPoint = 3;
  double* rloct, * weightt;
  GaussianQuadrature(qPoint,rloct,weightt);

  int itp = wParticlesPtr[prodWellMap[wInd]].start;
  for(int side=0;side<4;side++){
    if(wParticlesPtr[prodWellMap[wInd]].side[side]){
      for(int sec=0;sec<int_section_per_side;sec++){
	for(int qp=0;qp<int_quad_point;qp++){

	  double lntauMn,lntauVar;
	  double lntauMin,lntauMax,lntauFront,lPoint,rPoint;
	  int len_tmp = bTrackPtr->getParticle(itp)->getLength();
	  double tauMn = bTrackPtr->getParticle(itp)->getTrvlTimeAvg(len_tmp-1);
	  double tauVar = bTrackPtr->getParticle(itp)->getTrvlTimeVar(len_tmp-1);
	  int lenn = 1;
	  transform1(lenn,&tauMn,&tauVar,&lntauMn,&lntauVar);
	  lntauMin = lPoint = lntauMn - (sDRange*sqrt(lntauVar));
	  lntauMax = lntauMn + (sDRange*sqrt(lntauVar));
	  lntauFront = log(blSolnPtr->getFpwStar()*wTime);
	  while(true){
	    rPoint = lPoint+(intRange*sqrt(lntauVar));
	    if(rPoint>lntauMax){rPoint = lntauMax;}
	    if((lntauFront<rPoint)&&(lntauFront>lPoint)){rPoint = lntauFront;}
	    for(int iq = 0; iq < qPoint; iq++) {
	      double lntauVal = ((rPoint-lPoint)/2*rloct[iq])+((rPoint+lPoint)/2);
	      double tauVal = exp(lntauVal);
	      double fwVal = blSolnPtr->getFracFlow1(tauVal/wTime);
	      qwMn += pow(-1.,side)*lSection[side]/2*weightw[qp]*darcyV[side]*dz*
		      (rPoint-lPoint)/2*weightt[iq]*fwVal*pdfGauss(lntauVal,lntauMn,lntauVar);
	    }
	    if(rPoint==lntauMax){break;};
	    lPoint = rPoint;
	  }
	  itp++;
	}
      }
    }
  }
  delete[] rlocw;
  delete[] weightw;
  delete[] rloct;
  delete[] weightt;

  return qwMn;
}

double Production::CalcWaterProdRateCov(int wInd1, double wTime1, double qwMn1, int wInd2, double wTime2, double qwMn2){
  double qwCov = 0;
  Well* wellsPtr = domainPtr->getWellArr();
  double* rhoLntau = bTrackPtr->getLogTrvlTimeRho();

  //The current implementation allows only one-block well
  int wi1 = wellsPtr[wInd1].getWellIJK(0,0);
  int wj1 = wellsPtr[wInd1].getWellIJK(1,0);
  int wk1 = wellsPtr[wInd1].getWellIJK(2,0);
  int wi2 = wellsPtr[wInd2].getWellIJK(0,0);
  int wj2 = wellsPtr[wInd2].getWellIJK(1,0);
  int wk2 = wellsPtr[wInd2].getWellIJK(2,0);
  double dx1 = domainPtr->getBlock(wi1,wj1,wk1)->getDx();
  double dy1 = domainPtr->getBlock(wi1,wj1,wk1)->getDy();
  double dz1 = domainPtr->getBlock(wi1,wj1,wk1)->getDz();
  double poro1 = domainPtr->getPoro(wi1,wj1,wk1);
  double dx2 = domainPtr->getBlock(wi2,wj2,wk2)->getDx();
  double dy2 = domainPtr->getBlock(wi2,wj2,wk2)->getDy();
  double dz2 = domainPtr->getBlock(wi2,wj2,wk2)->getDz();
  double poro2 = domainPtr->getPoro(wi2,wj2,wk2);
  
  double darcyV1[4];
  darcyV1[0] = domainPtr->getBlock(wi1,wj1,wk1)->getVx1()*poro1;
  darcyV1[1] = domainPtr->getBlock(wi1,wj1,wk1)->getVx2()*poro1;
  darcyV1[2] = domainPtr->getBlock(wi1,wj1,wk1)->getVy1()*poro1;
  darcyV1[3] = domainPtr->getBlock(wi1,wj1,wk1)->getVy2()*poro1;
  double darcyV2[4];
  darcyV2[0] = domainPtr->getBlock(wi2,wj2,wk2)->getVx1()*poro2;
  darcyV2[1] = domainPtr->getBlock(wi2,wj2,wk2)->getVx2()*poro2;
  darcyV2[2] = domainPtr->getBlock(wi2,wj2,wk2)->getVy1()*poro2;
  darcyV2[3] = domainPtr->getBlock(wi2,wj2,wk2)->getVy2()*poro2;
  double lSection1[4] = {dy1/int_section_per_side,dy1/int_section_per_side,dx1/int_section_per_side,dx1/int_section_per_side};
  double lSection2[4] = {dy2/int_section_per_side,dy2/int_section_per_side,dx2/int_section_per_side,dx2/int_section_per_side};

  double* Cv1v1 = domainPtr->getCv1v1();
  double* Cv1v2 = domainPtr->getCv1v2();
  double* Cv2v2 = domainPtr->getCv2v2();

  double* rlocw, * weightw;
  GaussianQuadrature(int_quad_point,rlocw,weightw);
  double sDRange = 3.5, intRange = 0.5;
  int qPoint = 3;
  double* rloct, * weightt;
  GaussianQuadrature(qPoint,rloct,weightt);

  int itp1 = wParticlesPtr[prodWellMap[wInd1]].start;
  for(int side1=0;side1<4;side1++){
    if(wParticlesPtr[prodWellMap[wInd1]].side[side1]){	    
      for(int sec1=0;sec1<int_section_per_side;sec1++){      
	for(int qp1=0;qp1<int_quad_point;qp1++){

	  double lntauMn1,lntauVar1;
	  double lntauMin1,lntauMax1,lntauFront1,lPoint1,rPoint1;
	  int len_tmp1 = bTrackPtr->getParticle(itp1)->getLength();
	  double tauMn1 = bTrackPtr->getParticle(itp1)->getTrvlTimeAvg(len_tmp1-1);
	  double tauVar1 = bTrackPtr->getParticle(itp1)->getTrvlTimeVar(len_tmp1-1);
	  int lenn = 1;
	  transform1(lenn,&tauMn1,&tauVar1,&lntauMn1,&lntauVar1);
	  lntauMin1 = lPoint1 = lntauMn1 - (sDRange*sqrt(lntauVar1));
	  lntauMax1 = lntauMn1 + (sDRange*sqrt(lntauVar1));
	  lntauFront1 = log(blSolnPtr->getFpwStar()*wTime1);

	  while(true){
	    rPoint1 = lPoint1+(intRange*sqrt(lntauVar1));
	    if(rPoint1>lntauMax1){rPoint1 = lntauMax1;}
	    if((lntauFront1<rPoint1)&&(lntauFront1>lPoint1)){rPoint1 = lntauFront1;}
	    for(int iq1=0;iq1<qPoint;iq1++){
	      double lntauVal1 = ((rPoint1-lPoint1)/2*rloct[iq1])+((rPoint1+lPoint1)/2);
	      double tauVal1 = exp(lntauVal1);
	      double fwVal1 = blSolnPtr->getFracFlow1(tauVal1/wTime1);
	      
	      int itp2 = wParticlesPtr[prodWellMap[wInd2]].start;
	      for(int side2=0;side2<4;side2++){
		if(wParticlesPtr[prodWellMap[wInd2]].side[side2]){
		  
		  double Cvv;
		  if     ((side1==0)&&(side2==0)) Cvv = Cv1v1[domainPtr->getIndex(wi1-1,wj1,wk1)+domainPtr->getIndex(wi2-1,wj2,wk2)*gridPtr->getnNode()];
		  else if((side1==0)&&(side2==1)) Cvv = Cv1v1[domainPtr->getIndex(wi1-1,wj1,wk1)+domainPtr->getIndex(wi2,wj2,wk2)*gridPtr->getnNode()];
		  else if((side1==0)&&(side2==2)) Cvv = Cv1v2[domainPtr->getIndex(wi1-1,wj1,wk1)+domainPtr->getIndex(wi2,wj2-1,wk2)*gridPtr->getnNode()];
		  else if((side1==0)&&(side2==3)) Cvv = Cv1v2[domainPtr->getIndex(wi1-1,wj1,wk1)+domainPtr->getIndex(wi2,wj2,wk2)*gridPtr->getnNode()];
		  else if((side1==1)&&(side2==0)) Cvv = Cv1v1[domainPtr->getIndex(wi1,wj1,wk1)+domainPtr->getIndex(wi2-1,wj2,wk2)*gridPtr->getnNode()];
		  else if((side1==1)&&(side2==1)) Cvv = Cv1v1[domainPtr->getIndex(wi1,wj1,wk1)+domainPtr->getIndex(wi2,wj2,wk2)*gridPtr->getnNode()];
		  else if((side1==1)&&(side2==2)) Cvv = Cv1v2[domainPtr->getIndex(wi1,wj1,wk1)+domainPtr->getIndex(wi2,wj2-1,wk2)*gridPtr->getnNode()];
		  else if((side1==1)&&(side2==3)) Cvv = Cv1v2[domainPtr->getIndex(wi1,wj1,wk1)+domainPtr->getIndex(wi2,wj2,wk2)*gridPtr->getnNode()];
		  else if((side1==2)&&(side2==0)) Cvv = Cv1v2[domainPtr->getIndex(wi2-1,wj2,wk2)+domainPtr->getIndex(wi1,wj1-1,wk1)*gridPtr->getnNode()];
		  else if((side1==2)&&(side2==1)) Cvv = Cv1v2[domainPtr->getIndex(wi2,wj2,wk2)+domainPtr->getIndex(wi1,wj1-1,wk1)*gridPtr->getnNode()];
		  else if((side1==2)&&(side2==2)) Cvv = Cv2v2[domainPtr->getIndex(wi1,wj1-1,wk1)+domainPtr->getIndex(wi2,wj2-1,wk2)*gridPtr->getnNode()];
		  else if((side1==2)&&(side2==3)) Cvv = Cv2v2[domainPtr->getIndex(wi1,wj1-1,wk1)+domainPtr->getIndex(wi2,wj2,wk2)*gridPtr->getnNode()];
		  else if((side1==3)&&(side2==0)) Cvv = Cv1v2[domainPtr->getIndex(wi2-1,wj2,wk2)+domainPtr->getIndex(wi1,wj1,wk1)*gridPtr->getnNode()];
		  else if((side1==3)&&(side2==1)) Cvv = Cv1v2[domainPtr->getIndex(wi2,wj2,wk2)+domainPtr->getIndex(wi1,wj1,wk1)*gridPtr->getnNode()];
		  else if((side1==3)&&(side2==2)) Cvv = Cv2v2[domainPtr->getIndex(wi1,wj1,wk1)+domainPtr->getIndex(wi2,wj2-1,wk2)*gridPtr->getnNode()];
		  else if((side1==3)&&(side2==3)) Cvv = Cv2v2[domainPtr->getIndex(wi1,wj1,wk1)+domainPtr->getIndex(wi2,wj2,wk2)*gridPtr->getnNode()];
		  
		  
		  for(int sec2=0;sec2<int_section_per_side;sec2++){
		    for(int qp2=0;qp2<int_quad_point;qp2++){
		      
		      double lntauMn2,lntauVar2;
		      double lntauMin2,lntauMax2,lntauFront2,lPoint2,rPoint2;
		      int len_tmp2 = bTrackPtr->getParticle(itp2)->getLength();
		      double tauMn2 = bTrackPtr->getParticle(itp2)->getTrvlTimeAvg(len_tmp2-1);
		      double tauVar2 = bTrackPtr->getParticle(itp2)->getTrvlTimeVar(len_tmp2-1);
		      transform1(lenn,&tauMn2,&tauVar2,&lntauMn2,&lntauVar2);
		      lntauMin2 = lPoint2 = lntauMn2 - (sDRange*sqrt(lntauVar2));
		      lntauMax2 = lntauMn2 + (sDRange*sqrt(lntauVar2));
		      lntauFront2 = log(blSolnPtr->getFpwStar()*wTime2);
		      double rhoLntau12 = rhoLntau[itp1*nParticles+itp2];
		      if(rhoLntau12<-0.99999){rhoLntau12 = -0.99999;}//to avoid strong negative correlation		      
 
		      if(rhoLntau12>0.99999){ //assume calculating variance, not covariance
			qwCov += lSection1[side1]/2*weightw[qp1]*lSection2[side2]/2*weightw[qp2]*
			         (rPoint1-lPoint1)/2*weightt[iq1]*
			         (pow(darcyV1[side1]*dz1*fwVal1,2)+(Cvv*poro1*poro1*fwVal1*fwVal1*dz1*dz1))*
			         pdfGauss(lntauVal1,lntauMn1,lntauVar1);
		      }
		      else{
			while(true){
			  rPoint2 = lPoint2+(intRange*sqrt(lntauVar2));
			  if(rPoint2>lntauMax2){rPoint2 = lntauMax2;}
			  if((lntauFront2<rPoint2)&&(lntauFront2>lPoint2)){rPoint2 = lntauFront2;}
			  for(int iq2=0;iq2<qPoint;iq2++) {
			    double lntauVal2 = ((rPoint2-lPoint2)/2*rloct[iq2])+((rPoint2+lPoint2)/2);
			    double tauVal2 = exp(lntauVal2);
			    double fwVal2 = blSolnPtr->getFracFlow1(tauVal2/wTime2);
			    
			    qwCov += pow(-1.,side1)*pow(-1.,side2)*
			      lSection1[side1]/2*weightw[qp1]*lSection2[side2]/2*weightw[qp2]*
			      (rPoint1-lPoint1)/2*weightt[iq1]*(rPoint2-lPoint2)/2*weightt[iq2]*
			      ((darcyV1[side1]*dz1*fwVal1*darcyV2[side2]*dz2*fwVal2)+
			       (Cvv*poro1*poro2*fwVal1*fwVal2*dz1*dz2))*
			      pdfGauss2(lntauVal1,lntauMn1,lntauVar1,lntauVal2,lntauMn2,lntauVar2,rhoLntau12);
			  }
			  if(rPoint2==lntauMax2){break;};
			  lPoint2 = rPoint2;
			}
		      }
		      /*
		      if(qwCov!=qwCov){//Check NAN number
			cout<<rhoLntau12<<" "
			    <<side1<<" "<<side2<<" "
			    <<sec1<<" "<<sec2<<" "
			    <<qp1<<" "<<qp2<<" ";
			if(rhoLntau12==1)
			  cout<<pdfGauss(lntauVal1,lntauMn1,lntauVar1)<<endl;
			else
			  cout<<lntauMn1<<" "<<lntauMn2<<" "<<lntauVar1<<" "<<lntauVar2<<endl;
			exit(1);
		      }
		      */
		      itp2++;
		    }
		  }
		}
	      }
	    }
	    if(rPoint1==lntauMax1){break;};
	    lPoint1 = rPoint1;
	  }
	  itp1++;
	}
      }
    }
  }
  //if(wTime1>100) cout<<"qwCov: "<<qwCov<<" ";

  itp1 = wParticlesPtr[prodWellMap[wInd1]].start;
  for(int side1=0;side1<4;side1++){
    if(wParticlesPtr[prodWellMap[wInd1]].side[side1]){	    
      for(int sec1=0;sec1<int_section_per_side;sec1++){
	for(int qp1=0;qp1<int_quad_point;qp1++){
		  
	  double lntauMn,lntauVar;
	  double lntauMin,lntauMax,lntauFront,lPoint,rPoint;
	  int len_tmp = bTrackPtr->getParticle(itp1)->getLength();
	  double tauMn = bTrackPtr->getParticle(itp1)->getTrvlTimeAvg(len_tmp-1);
	  double tauVar = bTrackPtr->getParticle(itp1)->getTrvlTimeVar(len_tmp-1);
	  int lenn = 1;
	  transform1(lenn,&tauMn,&tauVar,&lntauMn,&lntauVar);
	  lntauMin = lPoint = lntauMn - (sDRange*sqrt(lntauVar));
	  lntauMax = lntauMn + (sDRange*sqrt(lntauVar));
	  lntauFront = log(blSolnPtr->getFpwStar()*wTime1);
	  while(true){
	    rPoint = lPoint+(intRange*sqrt(lntauVar));
	    if(rPoint>lntauMax){rPoint = lntauMax;}
	    if((lntauFront<rPoint)&&(lntauFront>lPoint)){rPoint = lntauFront;}
	    for(int iq = 0; iq < qPoint; iq++) {
	      double lntauVal = ((rPoint-lPoint)/2*rloct[iq])+((rPoint+lPoint)/2);
	      double tauVal = exp(lntauVal);
	      double fwVal;
	      fwVal = blSolnPtr->getFracFlow1(tauVal/wTime1);
	      qwCov -= pow(-1.,side1)*lSection1[side1]/2*weightw[qp1]*darcyV1[side1]*dz1*
		(rPoint-lPoint)/2*weightt[iq]*fwVal*qwMn2*pdfGauss(lntauVal,lntauMn,lntauVar);
	    }
	    if(rPoint==lntauMax){break;};
	    lPoint = rPoint;
	  }
	  itp1++;
	}
      }
    }
  }
  //if(wTime1>100) cout<<qwCov<<" ";

  int itp2 = wParticlesPtr[prodWellMap[wInd2]].start;
  for(int side2=0;side2<4;side2++){
    if(wParticlesPtr[prodWellMap[wInd2]].side[side2]){
      for(int sec2=0;sec2<int_section_per_side;sec2++){
	for(int qp2=0;qp2<int_quad_point;qp2++){
	  
	  
	  double lntauMn,lntauVar;
	  double lntauMin,lntauMax,lntauFront,lPoint,rPoint;
	  int len_tmp = bTrackPtr->getParticle(itp2)->getLength();
	  double tauMn = bTrackPtr->getParticle(itp2)->getTrvlTimeAvg(len_tmp-1);
	  double tauVar = bTrackPtr->getParticle(itp2)->getTrvlTimeVar(len_tmp-1);
	  int lenn = 1;
	  transform1(lenn,&tauMn,&tauVar,&lntauMn,&lntauVar);
	  lntauMin = lPoint = lntauMn - (sDRange*sqrt(lntauVar));
	  lntauMax = lntauMn + (sDRange*sqrt(lntauVar));
	  lntauFront = log(blSolnPtr->getFpwStar()*wTime2);
	  while(true){
	    rPoint = lPoint+(intRange*sqrt(lntauVar));
	    if(rPoint>lntauMax){rPoint = lntauMax;}
	    if((lntauFront<rPoint)&&(lntauFront>lPoint)){rPoint = lntauFront;}
	    for(int iq = 0; iq < qPoint; iq++) {
	      double lntauVal = ((rPoint-lPoint)/2*rloct[iq])+((rPoint+lPoint)/2);
	      double tauVal = exp(lntauVal);
	      double fwVal;
	      fwVal = blSolnPtr->getFracFlow1(tauVal/wTime2);
	      qwCov -= pow(-1.,side2)*lSection2[side2]/2*weightw[qp2]*darcyV2[side2]*dz2*
		(rPoint-lPoint)/2*weightt[iq]*fwVal*qwMn1*pdfGauss(lntauVal,lntauMn,lntauVar);
	    }
	    if(rPoint==lntauMax){break;};
	    lPoint = rPoint;
	  }
	  itp2++;
	}
      }
    }
  }
  //if(wTime1>100) cout<<qwCov<<" ";

  qwCov += qwMn1*qwMn2;
  //if(wTime1>100) cout<<qwCov<<endl;

  delete[] rlocw;
  delete[] weightw;
  delete[] rloct;
  delete[] weightt;

  return qwCov;
}

double Production::CalcTotal_WaterProdRateCov(int wTInd, int wwInd, double wwTime){
  double qTwCov = 0;
  Well* wellsPtr = domainPtr->getWellArr();
  double* rhoLntau = bTrackPtr->getLogTrvlTimeRho();

  //The current implementation allows only one-block well
  int wiT = wellsPtr[wTInd].getWellIJK(0,0);
  int wjT = wellsPtr[wTInd].getWellIJK(1,0);
  int wkT = wellsPtr[wTInd].getWellIJK(2,0);
  int wiw = wellsPtr[wwInd].getWellIJK(0,0);
  int wjw = wellsPtr[wwInd].getWellIJK(1,0);
  int wkw = wellsPtr[wwInd].getWellIJK(2,0);
  double dxT = domainPtr->getBlock(wiT,wjT,wkT)->getDx();
  double dyT = domainPtr->getBlock(wiT,wjT,wkT)->getDy();
  double dzT = domainPtr->getBlock(wiT,wjT,wkT)->getDz();
  double poroT = domainPtr->getPoro(wiT,wjT,wkT);
  double dxw = domainPtr->getBlock(wiw,wjw,wkw)->getDx();
  double dyw = domainPtr->getBlock(wiw,wjw,wkw)->getDy();
  double dzw = domainPtr->getBlock(wiw,wjw,wkw)->getDz();
  double porow = domainPtr->getPoro(wiw,wjw,wkw);
  
  double darcyVT[4];
  darcyVT[0] = domainPtr->getBlock(wiT,wjT,wkT)->getVx1()*poroT;
  darcyVT[1] = domainPtr->getBlock(wiT,wjT,wkT)->getVx2()*poroT;
  darcyVT[2] = domainPtr->getBlock(wiT,wjT,wkT)->getVy1()*poroT;
  darcyVT[3] = domainPtr->getBlock(wiT,wjT,wkT)->getVy2()*poroT;
  double darcyVw[4];
  darcyVw[0] = domainPtr->getBlock(wiw,wjw,wkw)->getVx1()*porow;
  darcyVw[1] = domainPtr->getBlock(wiw,wjw,wkw)->getVx2()*porow;
  darcyVw[2] = domainPtr->getBlock(wiw,wjw,wkw)->getVy1()*porow;
  darcyVw[3] = domainPtr->getBlock(wiw,wjw,wkw)->getVy2()*porow;
  double lSectionT[4] = {dyT,dyT,dxT,dxT};
  double lSectionw[4] = {dyw/int_section_per_side,dyw/int_section_per_side,dxw/int_section_per_side,dxw/int_section_per_side};

  double* Cv1v1 = domainPtr->getCv1v1();
  double* Cv1v2 = domainPtr->getCv1v2();
  double* Cv2v2 = domainPtr->getCv2v2();

  double* rlocw, * weightw;
  GaussianQuadrature(int_quad_point,rlocw,weightw);
  double sDRange = 3.5, intRange = 1;
  int qPoint = 3;
  double* rloct, * weightt;
  GaussianQuadrature(qPoint,rloct,weightt);

  for(int sideT=0;sideT<4;sideT++){
    if(wParticlesPtr[prodWellMap[wTInd]].side[sideT]){	    

      int itpw = wParticlesPtr[prodWellMap[wwInd]].start;
      for(int sidew=0;sidew<4;sidew++){
	if(wParticlesPtr[prodWellMap[wwInd]].side[sidew]){
		  
	  double Cvv;
	  if     ((sideT==0)&&(sidew==0)) Cvv = Cv1v1[domainPtr->getIndex(wiT-1,wjT,wkT)+domainPtr->getIndex(wiw-1,wjw,wkw)*gridPtr->getnNode()];
	  else if((sideT==0)&&(sidew==1)) Cvv = Cv1v1[domainPtr->getIndex(wiT-1,wjT,wkT)+domainPtr->getIndex(wiw,wjw,wkw)*gridPtr->getnNode()];
	  else if((sideT==0)&&(sidew==2)) Cvv = Cv1v2[domainPtr->getIndex(wiT-1,wjT,wkT)+domainPtr->getIndex(wiw,wjw-1,wkw)*gridPtr->getnNode()];
	  else if((sideT==0)&&(sidew==3)) Cvv = Cv1v2[domainPtr->getIndex(wiT-1,wjT,wkT)+domainPtr->getIndex(wiw,wjw,wkw)*gridPtr->getnNode()];
	  else if((sideT==1)&&(sidew==0)) Cvv = Cv1v1[domainPtr->getIndex(wiT,wjT,wkT)+domainPtr->getIndex(wiw-1,wjw,wkw)*gridPtr->getnNode()];
	  else if((sideT==1)&&(sidew==1)) Cvv = Cv1v1[domainPtr->getIndex(wiT,wjT,wkT)+domainPtr->getIndex(wiw,wjw,wkw)*gridPtr->getnNode()];
	  else if((sideT==1)&&(sidew==2)) Cvv = Cv1v2[domainPtr->getIndex(wiT,wjT,wkT)+domainPtr->getIndex(wiw,wjw-1,wkw)*gridPtr->getnNode()];
	  else if((sideT==1)&&(sidew==3)) Cvv = Cv1v2[domainPtr->getIndex(wiT,wjT,wkT)+domainPtr->getIndex(wiw,wjw,wkw)*gridPtr->getnNode()];
	  else if((sideT==2)&&(sidew==0)) Cvv = Cv1v2[domainPtr->getIndex(wiw-1,wjw,wkw)+domainPtr->getIndex(wiT,wjT-1,wkT)*gridPtr->getnNode()];
	  else if((sideT==2)&&(sidew==1)) Cvv = Cv1v2[domainPtr->getIndex(wiw,wjw,wkw)+domainPtr->getIndex(wiT,wjT-1,wkT)*gridPtr->getnNode()];
	  else if((sideT==2)&&(sidew==2)) Cvv = Cv2v2[domainPtr->getIndex(wiT,wjT-1,wkT)+domainPtr->getIndex(wiw,wjw-1,wkw)*gridPtr->getnNode()];
	  else if((sideT==2)&&(sidew==3)) Cvv = Cv2v2[domainPtr->getIndex(wiT,wjT-1,wkT)+domainPtr->getIndex(wiw,wjw,wkw)*gridPtr->getnNode()];
	  else if((sideT==3)&&(sidew==0)) Cvv = Cv1v2[domainPtr->getIndex(wiw-1,wjw,wkw)+domainPtr->getIndex(wiT,wjT,wkT)*gridPtr->getnNode()];
	  else if((sideT==3)&&(sidew==1)) Cvv = Cv1v2[domainPtr->getIndex(wiw,wjw,wkw)+domainPtr->getIndex(wiT,wjT,wkT)*gridPtr->getnNode()];
	  else if((sideT==3)&&(sidew==2)) Cvv = Cv2v2[domainPtr->getIndex(wiT,wjT,wkT)+domainPtr->getIndex(wiw,wjw-1,wkw)*gridPtr->getnNode()];
	  else if((sideT==3)&&(sidew==3)) Cvv = Cv2v2[domainPtr->getIndex(wiT,wjT,wkT)+domainPtr->getIndex(wiw,wjw,wkw)*gridPtr->getnNode()];
		  
	  for(int secw=0;secw<int_section_per_side;secw++){
	    for(int qpw=0;qpw<int_quad_point;qpw++){
		      
	      double lntauMn,lntauVar;
	      double lntauMin,lntauMax,lntauFront,lPoint,rPoint;
	      int len_tmp = bTrackPtr->getParticle(itpw)->getLength();
	      double tauMn = bTrackPtr->getParticle(itpw)->getTrvlTimeAvg(len_tmp-1);
	      double tauVar = bTrackPtr->getParticle(itpw)->getTrvlTimeVar(len_tmp-1);
	      int lenn = 1;
	      transform1(lenn,&tauMn,&tauVar,&lntauMn,&lntauVar);
	      lntauMin = lPoint = lntauMn - (sDRange*sqrt(lntauVar));
	      lntauMax = lntauMn + (sDRange*sqrt(lntauVar));
	      lntauFront = log(blSolnPtr->getFpwStar()*wwTime);
	      
	      while(true){
		rPoint = lPoint+(intRange*sqrt(lntauVar));
		if(rPoint>lntauMax){rPoint = lntauMax;}
		if((lntauFront<rPoint)&&(lntauFront>lPoint)){rPoint = lntauFront;}
		for(int iq=0;iq<qPoint;iq++) {
		  double lntauVal = ((rPoint-lPoint)/2*rloct[iq])+((rPoint+lPoint)/2);
		  double tauVal = exp(lntauVal);
		  double fwVal = blSolnPtr->getFracFlow1(tauVal/wwTime);
		  
		  qTwCov += pow(-1.,sideT)*pow(-1.,sidew)*
		           lSectionT[sideT]*lSectionw[sidew]/2*weightw[qpw]*
		           (rPoint-lPoint)/2*weightt[iq]*
		           (Cvv*poroT*porow*fwVal*dzT*dzw)*
		           pdfGauss(lntauVal,lntauMn,lntauVar);
		}
		if(rPoint==lntauMax){break;};
		lPoint = rPoint;
	      }
	      itpw++;
	    }
	  }
	}
      }
    }
  }

  delete[] rlocw;
  delete[] weightw;
  delete[] rloct;
  delete[] weightt;

  return qTwCov;
}

double Production::CalcCrossCovQwY(int wInd, double wTime, int yInd){
  double qwYCov = 0;
  Well* wellsPtr = domainPtr->getWellArr();
  Perm* permPtr = flowPtr->getPermPtr();
  if(covYLnTau[prodWellMap[wInd]]==NULL) CalcCovYLnTau(wInd);
  
  SwCondiInverse* swCondPtr = new SwCondiInverse(domainPtr,flowPtr);
  swCondPtr->RefreshPressureGradient();

  //The current implementation allows only one-block well
  int wi = wellsPtr[wInd].getWellIJK(0,0);
  int wj = wellsPtr[wInd].getWellIJK(1,0);
  int wk = wellsPtr[wInd].getWellIJK(2,0);
  double dx = domainPtr->getBlock(wi,wj,wk)->getDx();
  double dy = domainPtr->getBlock(wi,wj,wk)->getDy();
  double dz = domainPtr->getBlock(wi,wj,wk)->getDz();
  double poro = domainPtr->getPoro(wi,wj,wk);
  
  double darcyV[4];
  darcyV[0] = domainPtr->getBlock(wi,wj,wk)->getVx1()*poro;
  darcyV[1] = domainPtr->getBlock(wi,wj,wk)->getVx2()*poro;
  darcyV[2] = domainPtr->getBlock(wi,wj,wk)->getVy1()*poro;
  darcyV[3] = domainPtr->getBlock(wi,wj,wk)->getVy2()*poro;
  double lSection[4] = {dy/int_section_per_side,dy/int_section_per_side,dx/int_section_per_side,dx/int_section_per_side};

  double* rlocw, * weightw;
  GaussianQuadrature(int_quad_point,rlocw,weightw);
  double sDRange = 3.5, intRange = 0.5;
  int qPoint = 3;
  double* rlocty, * weightty;
  GaussianQuadrature(qPoint,rlocty,weightty);

  double yavg = permPtr->getYAvg(yInd);
  double yvar = permPtr->getYVar(yInd);
  double ymin = yavg - (sDRange*sqrt(yvar));
  double ymax = yavg + (sDRange*sqrt(yvar));

  double cvy[4];
  cvy[0] = swCondPtr->CalcCovYVx(wi-1,wj,wk,wi,wj,wk,yInd);
  cvy[1] = swCondPtr->CalcCovYVx(wi,wj,wk,wi+1,wj,wk,yInd);
  cvy[2] = swCondPtr->CalcCovYVy(wi,wj-1,wk,wi,wj,wk,yInd);
  cvy[3] = swCondPtr->CalcCovYVy(wi,wj,wk,wi,wj+1,wk,yInd);

  int itp = wParticlesPtr[prodWellMap[wInd]].start;
  int itpstart = itp;
  for(int side=0;side<4;side++){
    if(wParticlesPtr[prodWellMap[wInd]].side[side]){
      double intfw = 0;
      for(int sec=0;sec<int_section_per_side;sec++){
	for(int qp=0;qp<int_quad_point;qp++){

	  double lntauMn,lntauVar;
	  double lntauMin,lntauMax,lntauFront,lPoint,rPoint;
	  int len_tmp = bTrackPtr->getParticle(itp)->getLength();
	  double tauMn = bTrackPtr->getParticle(itp)->getTrvlTimeAvg(len_tmp-1);
	  double tauVar = bTrackPtr->getParticle(itp)->getTrvlTimeVar(len_tmp-1);
	  int lenn = 1;
	  transform1(lenn,&tauMn,&tauVar,&lntauMn,&lntauVar);
	  double rhoYLnTau = covYLnTau[prodWellMap[wInd]][(itp-itpstart)*gridPtr->getNumPermNode()+yInd]/sqrt(yvar*lntauVar);

	  lntauMin = lPoint = lntauMn - (sDRange*sqrt(lntauVar));
	  lntauMax = lntauMn + (sDRange*sqrt(lntauVar));
	  lntauFront = log(blSolnPtr->getFpwStar()*wTime);
	  while(true){
	    rPoint = lPoint+(intRange*sqrt(lntauVar));
	    if(rPoint>lntauMax){rPoint = lntauMax;}
	    if((lntauFront<rPoint)&&(lntauFront>lPoint)){rPoint = lntauFront;}
	    for(int iqt = 0; iqt < qPoint; iqt++) {
	      double lntauVal = ((rPoint-lPoint)/2*rlocty[iqt])+((rPoint+lPoint)/2);
	      double tauVal = exp(lntauVal);
	      double fwVal = blSolnPtr->getFracFlow1(tauVal/wTime);

	      intfw += lSection[side]/2*weightw[qp]*(rPoint-lPoint)/2*weightty[iqt]*
		       fwVal*pdfGauss(lntauVal,lntauMn,lntauVar);

	      double lPointY,rPointY;
	      lPointY = ymin;
	      while(true){
	        rPointY = lPointY+(intRange*sqrt(yvar));
		if(rPointY>ymax) rPointY = ymax;
		for(int iqy=0;iqy<qPoint;iqy++){
		  double yval = ((rPointY-lPointY)/2*rlocty[iqy])+((rPointY+lPointY)/2);
		  
		  qwYCov += pow(-1.,side)*lSection[side]/2*weightw[qp]*
		            (rPoint-lPoint)/2*weightty[iqt]*(rPointY-lPointY)/2*weightty[iqy]*
			    darcyV[side]*dz*fwVal*(yval-yavg)*
			    pdfGauss2(lntauVal,lntauMn,lntauVar,yval,yavg,yvar,rhoYLnTau);
		}
		if(rPointY==ymax){break;};
		lPointY = rPointY;
	      }
	    }
	    if(rPoint==lntauMax){break;};
	    lPoint = rPoint;
	  }
	  itp++;
	}
      }
      qwYCov += pow(-1.,side)*cvy[side]*poro*dz*intfw;
    }
  }
  delete[] rlocw;
  delete[] weightw;
  delete[] rlocty;
  delete[] weightty;

  delete swCondPtr;

  return qwYCov;
}

