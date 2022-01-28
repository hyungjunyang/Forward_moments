#include "PCondiInverse.h"

PCondiInverse::PCondiInverse(Grid* grid, Perm *perm, P0Eqn* P0ee, CYPEqn *CYPee, 
                             CPPEqn *CPPee, bool debug_) 
 : gridPtr(grid), permPtr(perm), P0e( P0ee ), CYPe (CYPee), CPPe (CPPee),
   debug( debug_ )
{  
   if(debug) {
      cout<< endl 
          << " PCondiInverse::PCondiInverse()"<<endl;
   }
   
   p0  = P0e->getP0();

   readData();
   addArrays();
   addKrigWeighObj();
}

PCondiInverse::~PCondiInverse(){
   if(debug) 
      cout<<"~PCondiInverse::PCondiInverse()"<<endl;

   delKrigWeighObj();
   delArrays();
   delMem_of_readData();
}

void PCondiInverse::addMem_of_readData() {
   p_avg_Meas = new double[num_p_Meas];
   p_std_Meas = new double[num_p_Meas];
     i_p_Meas = new int   [num_p_Meas];
     j_p_Meas = new int   [num_p_Meas];
     k_p_Meas = new int   [num_p_Meas];
   ijk_p_Meas = new int   [num_p_Meas];
}

void PCondiInverse::delMem_of_readData() {
   delete [] p_avg_Meas;
   delete [] p_std_Meas;
   delete []   i_p_Meas;
   delete []   j_p_Meas;
   delete []   k_p_Meas;
   delete [] ijk_p_Meas;
}

void PCondiInverse::addArrays() {
   p_avg_UnCd = new double [ gridPtr->getnNode() ];
   p_std_UnCd = new double [ gridPtr->getnNode() ];
   p_avg_Cond = new double [ gridPtr->getnNode() ];
   p_std_Cond = new double [ gridPtr->getnNode() ];
   p_avg_Filt = new double [ gridPtr->getnNode() ];
   p_std_Filt = new double [ gridPtr->getnNode() ];
   
   Y_avg_UnCd = new double [ gridPtr->getNumPermNode() ];
   Y_std_UnCd = new double [ gridPtr->getNumPermNode() ];
   Y_avg_Cond = new double [ gridPtr->getNumPermNode() ];
   Y_std_Cond = new double [ gridPtr->getNumPermNode() ];   
}

void PCondiInverse::delArrays() {
   delete [] p_avg_UnCd;
   delete [] p_std_UnCd;
   delete [] p_avg_Cond;
   delete [] p_std_Cond;
   delete [] p_avg_Filt;
   delete [] p_std_Filt;   

   delete [] Y_avg_UnCd;
   delete [] Y_std_UnCd;
   delete [] Y_avg_Cond;
   delete [] Y_std_Cond;
}

void PCondiInverse::readData() {
   ifstream in_Pcondi;
   in_Pcondi.open("press_condi.in",ios::in);
   if(in_Pcondi.bad() ) {
      cerr << " Can not open file " <<  "press_condi.in"  << endl;
      exit(8);
   } else {
      if(debug) 
         cout << "File press_condi.in was opened successfully!" <<endl;
   }
   
   //1) number of pressure measurements
   Junk(in_Pcondi); in_Pcondi >> num_p_Meas;
   if(debug) cout<<"num_p_Meas = " << num_p_Meas <<endl;

   //2) allocate memory
   addMem_of_readData();

   //3) Continue reading measurements
   for(int m = 0; m < num_p_Meas; ++m) {
       Junk(in_Pcondi); 
       in_Pcondi >> i_p_Meas[m] >> j_p_Meas[m] >> k_p_Meas[m]
                 >> p_avg_Meas[m];
       --i_p_Meas[m];
       --j_p_Meas[m];
       --k_p_Meas[m];
       ijk_p_Meas[m] = gridPtr->getIndex(i_p_Meas[m], j_p_Meas[m], k_p_Meas[m]);
       p_std_Meas[m] = 0.;
       if(debug) cout <<   i_p_Meas[m]<<' '
                      <<   j_p_Meas[m]<<' '
                      <<   k_p_Meas[m]<<' '
                      << ijk_p_Meas[m]<<' '
                      << p_avg_Meas[m]<<' '
                      << p_std_Meas[m]<< endl;
   }
   Junk(in_Pcondi); 
        in_Pcondi >> max_iter >> eps1 >> relax;
	//if(debug) 
           cout << "Iteration Control : "
		<< max_iter <<' '
	        << eps1  << ' '
	        << relax << ' ' << endl;
}

// ===== Section 3 Solve() =====
void PCondiInverse::solve(int nWells, const Well *wellArr, const Control &c, 
                           Solver &s, int time_step ) {
   if(debug) cout << "PCondiInverse::solve()" << endl;

   //1) set both Unconditonal Y and P moments in this object!;
   //   and write UnCd Pressure Standard Dev.
   if(debug) cout << "step 1" << endl;
   CPPe->calcUnCdPStd();
   setUnCondMoments();
   output_p("UnCondP", p_avg_UnCd, p_std_UnCd);
   permPtr->wrtPermMomnt("UnCdPerm");
   
   //2) write out the measurements data 
   if(debug) cout << "step 2" << endl;
   output_p_Meas("MeaSurP");
   
   //3) cokrigging
   if(debug) cout << "step 3" << endl; 

   double * mu   = new double[ gridPtr->getNumPermNode() * num_p_Meas];
   double * aa   = new double[ num_p_Meas * num_p_Meas ];
   double * aaInv= new double[ num_p_Meas * num_p_Meas ];
   calcMu(aa, aaInv, mu);
   // Set or Update Y_Avg
   double value;
   for(int k2 = 0; k2 < gridPtr->getPermNz(); ++k2) {
       for(int j2 = 0; j2 < gridPtr->getPermNy(); ++j2) {
           for(int i2 = 0; i2 < gridPtr->getPermNx(); ++i2) {
               int ijk2 = gridPtr->getPermIndex(i2, j2, k2);
               value = permPtr->getYAvg(ijk2);
               for(int msu_j = 0; msu_j < num_p_Meas; ++msu_j) {
                   value += mu[ msu_j + ijk2 * num_p_Meas ] 
                          * ( p_avg_Meas[msu_j] - p_avg_UnCd[ ijk_p_Meas[msu_j] ] );
               }
               permPtr->setYAvg(ijk2, value);
            }
       }
   }
   // Set C_YY
   for(int ijk2 = 0; ijk2 < gridPtr->getNumPermNode(); ++ijk2) {
       for(int ijk1 = 0; ijk1 < gridPtr->getNumPermNode(); ++ijk1) {
           value = 0.;
           for(int msu_j = 0; msu_j < num_p_Meas; ++msu_j) {
               value += mu[ msu_j + ijk2 * num_p_Meas ]
                        * CYPe->getCYP(ijk_p_Meas[msu_j], ijk1);
           }
           permPtr->setCYY_Cond(ijk1, ijk2, value);
       }	   
   }
   delete [] aa;
   delete [] aaInv;
   delete [] mu;
   // ForWard Calculations
    P0e->solve(nWells, wellArr, c, s, 0);
   CYPe->solve(nWells, wellArr, c, s, 0);
   CPPe->solve(nWells, wellArr, c, s, 0);
   CPPe->calcCondPStd();
   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_Cond[ i ] = p0[i];
       p_std_Cond[ i ] = CPPe->getCondPStd( i );
   }
   output_p("ConditP", p_avg_Cond, p_std_Cond);
   permPtr->wrtPermMomnt("CondPerm");

   //4) calculate CPP conditional to Pressure measurements, we
   //   have two algorithm(), one is linear interpretation,
   //   the other is the fixed pressure
   if(debug) cout << "step 4" << endl;
   int algorithm = 2;
   if( algorithm == 1) {
       condCPPbyGaussDist();
   } else {
       condCPPbySolveEqn(nWells, wellArr, c, s);
   }
   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_Filt[i] = p0[i];
       p_std_Filt[i] = CPPe->getCondPStd( i );
   }
   output_p("TarGetP", p_avg_Filt, p_std_Filt);
}

// ===== Section 4 SolveIter() =====
void PCondiInverse::solveIter(int nWells, const Well *wellArr, const Control &c, 
                              Solver &s, int time_step ) {
   debug = false;
   if(debug) cout << "PCondiInverse::solve()" << endl;
   
   //1) set both Unconditonal Y and P moments in this object!;
   //   and write UnCd Pressure Standard Dev.
   if(debug) cout << "step 1" << endl;
   CPPe->calcUnCdPStd();
   setUnCondMoments();
   output_p("UnCondP", p_avg_UnCd, p_std_UnCd);
   permPtr->wrtPermMomnt("UnCdPerm");
   
   //2) write out the measurements data 
   if(debug) cout << "step 2" << endl;
   output_p_Meas("MeaSurP");
   
   //3) Solve with Iteration Scheme
   if(debug) cout << "step 3" << endl;
   double * mu   = new double[ gridPtr->getNumPermNode() * num_p_Meas];
   double * aa   = new double[ num_p_Meas * num_p_Meas ];
   double * aaInv= new double[ num_p_Meas * num_p_Meas ];
  
   //3.0) Set the fields related to Iteration Scheme
   setIterFields();
   
   //3.1) Copy the unconditonal Y_Avg to a Tmp Array in case.
   permPtr->Copy_PermYAvgToTmp();

   double old_error, new_error;
   old_error = getDiff(p_avg_UnCd, error_old);
   //calcMu(aa, aaInv, mu);
   for(int iter = 0; iter < max_iter; ++iter) {
       if(fabs(old_error) < eps1) {
          //4.4) Update C_YY
		  double value;
          for(int ijk2 = 0; ijk2 < gridPtr->getNumPermNode(); ++ijk2) {
	      for(int ijk1 = 0; ijk1 < gridPtr->getNumPermNode(); ++ijk1) {
                  value = 0.;
                  for(int msu_j = 0; msu_j < num_p_Meas; ++msu_j) {
                      value += mu[ msu_j + ijk2 * num_p_Meas ]
                            * CYPe->getCYP(ijk_p_Meas[msu_j], ijk1);
                  }
				  permPtr->setCYY_Cond(ijk1, ijk2, value);
	      }	   
          }
	       
          //4.5) ForWard Calculations
          CYPe->solve(nWells, wellArr, c, s, 0);
          CPPe->solve(nWells, wellArr, c, s, 0);
          CPPe->calcCondPStd();
          for(int i = 0; i < gridPtr->getnNode(); ++i ) {
              p_std_Cond[ i ] = CPPe->getCondPStd( i );
          }
	      break;
       }
       cout << "Iter = " << iter <<" Iteration Results: " << endl;

       if(iter == 0) {
          for(int i = 0; i < num_p_Meas; ++i) {
              error_crt[i] = p_avg_Meas[i] - p_avg_UnCd[ ijk_p_Meas[i] ];
	      error_crt[i] *= relax;
              /*cout << p_avg_Meas[i] << ' '
	           << p_avg_UnCd[ ijk_p_Meas[i] ] << ' '
	           << error_crt[i] << endl;*/
          }
       } else{
          for(int i = 0; i < num_p_Meas; ++i) {
              error_crt[i] = p_avg_Meas[i] - p_avg_Cond[ ijk_p_Meas[i] ];
              /*cout << p_avg_Meas[i] << ' '
	           << p_avg_Cond[ ijk_p_Meas[i] ] << ' '
	           << error_crt[i] << endl;*/
			  error_crt[i] *= relax;
		  }
       }
       cout << "old error  = " << old_error << endl;
       cout << endl;
       
       //4.1) using C_PP and C_YP to calculate mu only!
       calcMu(aa, aaInv, mu);
          
       //4.2) using mu to update the mean Y and Mean K
       double value;
       for(int i = 0; i < gridPtr->getNumPermNode(); ++i) {
		   value = permPtr->getYAvg(i);
           for(int msu_j = 0; msu_j < num_p_Meas; ++msu_j) {
               value += mu[ msu_j + i * num_p_Meas ] * error_crt [ msu_j]; 
           }
           permPtr->setYAvg(i, value);           	   
       }
       //4.3) Solve for P0       
       P0e->solve(nWells, wellArr, c, s, 0);
       for(int i = 0; i < gridPtr->getnNode(); ++i ) {
           p_avg_Cond[ i ] = p0[i];
       }
       new_error = getDiff(p_avg_Cond, error_new);
       cout << "new error  = " << new_error << endl;
       //fabs( old_error - new_error ) < 0.00001 ||
       if(iter == max_iter - 1) {
	  //if( (iter+1) % 5 == 0){     
	  cout << "do cyy update at iter " << iter << endl;
       
          //4.4) Update C_YY
          for(int ijk2 = 0; ijk2 < gridPtr->getNumPermNode(); ++ijk2) {
	      for(int ijk1 = 0; ijk1 < gridPtr->getNumPermNode(); ++ijk1) {
                  value = 0.;
                  for(int msu_j = 0; msu_j < num_p_Meas; ++msu_j) {
                      value += mu[ msu_j + ijk2 * num_p_Meas ]
                            * CYPe->getCYP(ijk_p_Meas[msu_j], ijk1);
                  }
	          permPtr->setCYY_Cond(ijk1, ijk2, value);
	      }	   
          }
	       
          //4.5) ForWard Calculations
          CYPe->solve(nWells, wellArr, c, s, 0);
          CPPe->solve(nWells, wellArr, c, s, 0);
          CPPe->calcCondPStd();
          for(int i = 0; i < gridPtr->getnNode(); ++i ) {
              p_std_Cond[ i ] = CPPe->getCondPStd( i );
          }
	  //calcMu(aa, aaInv, mu);
       }
       old_error = new_error;
       /*
       for(int i = 0; i < num_p_Meas; ++i ) {
	   error_crt[i] = p_avg_Meas[i] - p_avg_Cond[ ijk_p_Meas[i] ];
           cout << p_avg_Meas[i] << ' '
	        << p_avg_Cond[ ijk_p_Meas[i] ] << ' '
	        << error_crt[i] << endl;
       }
       cout << endl;
       */
       
      
       /*
       for(int i = 0; i < num_p_Meas; ++i) 
           error_old[i] = error_crt[i];
       */
       /* 
       if(new_error > old_error) {
          cout<<"Iteration starts deterioration!!" << endl;
	  break;
       } else {
          old_error = new_error;
	  for(int i = 0; i < num_p_Meas; ++i) 
              error_old[i] = error_crt[i];
       }*/
   }
   delIterFields();
   delete [] aa;
   delete [] aaInv;
   delete [] mu;
  
   output_p("ConditP", p_avg_Cond, p_std_Cond);
   permPtr->wrtPermMomnt("CondPerm");
   
   //5) calculate CPP conditional to Pressure measurements, we
   //   have two algorithm(), one is linear interpretation,
   //   the other is the fixed pressure
   if(debug) cout << "step 4" << endl;
   int algorithm = 2;
   if( algorithm == 1) {
       condCPPbyGaussDist();
   } else {
       //permPtr->Copy_TmpToPermYAvg();	   
       //condCPPbySolveEqn(nWells, wellArr, c, s);//[Pipatl] comment out, not to impose measurement to the problem
   }
   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_Filt[i] = p0[i];
       p_std_Filt[i] = CPPe->getCondPStd( i );
   }
   output_p("TarGetP", p_avg_Filt, p_std_Filt);
}
//4.2 Set or Update Y_Avg
//permPtr->Copy_TmpToPermYAvg();

// using C_PP and C_YP
void PCondiInverse::calcMu(double *aa, double *aaInv, double *mu) {
   for(int msu_k = 0; msu_k < num_p_Meas; ++msu_k) {
       for(int msu_j = 0; msu_j < num_p_Meas; ++msu_j) {
           aa[msu_j + msu_k * num_p_Meas] = CPPe->getCPP(ijk_p_Meas[msu_j], ijk_p_Meas[msu_k]);
       }
   }
   krigPtr->calCovInverse(num_p_Meas, aa, aaInv); 

   if(debug) {
      cout << "Inverse Checking !!" << endl;
      double sum;
      cout << setiosflags(ios::fixed | ios::showpoint) 
           << setprecision(5);

      for(int msu_k = 0; msu_k < num_p_Meas; ++msu_k) {
          for(int msu_j = 0; msu_j < num_p_Meas; ++msu_j) {
              sum = 0;
              for(int k = 0; k < num_p_Meas; ++k) { 
                  sum += aa[msu_j + k * num_p_Meas] * aaInv[k + msu_k * num_p_Meas];
              }
              cout << sum << ' ';
          }
          cout << endl;
      }
   }
   // calculate mu - the weighting array
   double sum;
   for(int ijk = 0; ijk < gridPtr->getNumPermNode(); ++ijk) {
       for(int msu_k = 0; msu_k < num_p_Meas; ++msu_k) {
           sum = 0.;
           for(int msu_j = 0; msu_j < num_p_Meas; ++msu_j) { 
               sum += aaInv[msu_j + msu_k * num_p_Meas] * CYPe->getCYP(ijk_p_Meas[msu_j], ijk);
           }
           mu[msu_k + ijk * num_p_Meas] = sum;
       }
   }
}

void PCondiInverse::setIterFields(){
   //max_iter = 10;
   //eps1  = 0.01;
   //relax = 0.1;
   error_old = new double[num_p_Meas];
   error_crt = new double[num_p_Meas];
   error_new = new double[num_p_Meas];
   for(int i = 0; i < num_p_Meas; ++i) {
       error_old[i] = 0.;
       error_crt[i] = 0.;
       error_new[i] = 0.;
   }
}

void PCondiInverse::delIterFields(){
   delete [] error_old;
   delete [] error_crt;
   delete [] error_new;
}

double PCondiInverse::getDiff(double *p_avg_Calc, double *error_Calc) {
  double sum = 0.;
  for(int i = 0; i < num_p_Meas; ++i) {
      error_Calc[i] =  p_avg_Meas[i] - p_avg_Calc[ ijk_p_Meas[i] ];
      //cout << p_avg_Meas[i] << ' '
      //     << p_avg_Calc[ ijk_p_Meas[i] ] << ' '
      //     << error_Calc[i] << endl;
      sum += fabs(error_Calc[i]);
  }
  //cout << endl;
  sum = sum/num_p_Meas;
  return sum;
}

// ===== setUnCondMoments() =====
void PCondiInverse::setUnCondMoments() {
   for(int i = 0; i < gridPtr->getNumPermNode(); ++i ) {
       Y_avg_UnCd[i] =       permPtr->getYAvg( i );
       Y_std_UnCd[i] = sqrt( permPtr->getYVar( i ) ); 
   }
   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_UnCd[i] =  p0[i];
       p_std_UnCd[i] =  CPPe->getUnCdPStd( i );
   }
}

// ===== condCPPbySolveEqn()  =====
void PCondiInverse::condCPPbySolveEqn(int nWells, const Well *w, 
                                      const Control &c, Solver &s) {
   P0e->solve(nWells, w, c, s, 0, num_p_Meas, 
              i_p_Meas, j_p_Meas, k_p_Meas, p_avg_Meas);
   CYPe->solve(nWells, w, c, s, 0);
   CPPe->solve(nWells, w, c, s, 0);
   CPPe->calcCondPStd();
}

// ===== condCPPbyGaussDist()  =====
void PCondiInverse::condCPPbyGaussDist() {
   krigPtr->condGaussDist(gridPtr->getnNode(), num_p_Meas, 
                          ijk_p_Meas, CPPe->getCPP(),p_std_Filt);
}

void PCondiInverse::addKrigWeighObj() {
   krigPtr = new KrigWeigh(gridPtr, permPtr, num_p_Meas);
}

void PCondiInverse::delKrigWeighObj() {
   delete krigPtr;
}

////////////////////////////////////////////////////////////////////////
//////////////////////// OUTPUT ZONE ////////////////////////////////////
//===== output_p() =====
void PCondiInverse::output_p(char *file, double *avg, double *std){
   char * ending1 = "_1d.out";
   char * ending2 = "_2d.out";
   char * ending3 = "_2d.in";
   
   int length  = strlen(file);
   int length1 = strlen(ending1);
   int length2 = strlen(ending2);
   int length3 = strlen(ending2);
   
   char *file_1d   = new char [length + length1 + 2];   
   char *file_2d   = new char [length + length2 + 2];
   char *file_2din = new char [length + length3 + 2];
   
   strcpy(file_1d, file);
   strcpy(file_2d, file);
   strcpy(file_2din, file);
   
   strcat(file_1d,ending1);
   strcat(file_2d,ending2); 
   strcat(file_2din,ending3);
   
   ofstream os1(file_1d, ios::out);
   ofstream os2(file_2d, ios::out);
   ofstream os3(file_2din, ios::out);
   
   os1 << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);
   os2 << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);    
   os3 << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);
   for(int k = 0; k < gridPtr->getNz(); k++) {
       for(int j = 0; j < gridPtr->getNy(); j++) {
           int i = j;
           int ijk = gridPtr->getIndex(i, j, k);
           double dist = sqrt(  gridPtr->getX(i) * gridPtr->getX(i) 
                              + gridPtr->getY(i) * gridPtr->getY(i) );
           os1 << dist << ' '<< avg[ijk] << ' ' << std[ijk] << endl;
               
           for(int i = 0; i < gridPtr->getNx(); i++) {
               int ijk = gridPtr->getIndex(i, j, k);
               os2 << gridPtr->getX(i) << ' ' << gridPtr->getY(j) << ' '
                   << avg[ijk] << ' ' << std[ijk] << endl;
               os3 << i + 1 << ' ' << j + 1 << ' ' << k + 1<< ' '
                   << avg[ijk] << endl; 
           }
           os2<<endl;
       }
   }
   os1<<endl;
   os2<<endl;
   os3<<endl;
   
   os1.close();
   os2.close();
   os3.close();
   
   delete [] file_1d;
   delete [] file_2d;
   delete [] file_2din;
}

//===== output_p_Meas() =====
void PCondiInverse::output_p_Meas(char *file) {
   char * ending1 = "_1d.out";
   char * ending2 = "_2d.out";
   int length  = strlen(file);
   int length1 = strlen(ending1);
   int length2 = strlen(ending2);
   char *file_1d   = new char [length + length1 + 2];   
   char *file_2d   = new char [length + length2 + 2];
   
   strcpy(file_1d, file);
   strcpy(file_2d, file);
   
   strcat(file_1d,   ending1);
   strcat(file_2d,   ending2); 
   
   ofstream os1(file_1d, ios::out);
   ofstream os2(file_2d, ios::out);
   os1 << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);
   os2 << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);    
   for(int i = 0; i < num_p_Meas; ++i) {
       os2 << gridPtr->getX( i_p_Meas[i] ) << ' ' 
           << gridPtr->getY( j_p_Meas[i] ) << ' '
           << p_avg_Meas[i] << ' ' 
           << p_std_Meas[i] << endl;
 
       if( i_p_Meas[i] == j_p_Meas[i] ) {
           double dist = sqrt(  gridPtr->getX( i_p_Meas[i] ) 
                              * gridPtr->getX( i_p_Meas[i] )
                              + gridPtr->getY( j_p_Meas[i] )
                              * gridPtr->getY( j_p_Meas[i] )
                             );
           os1 << dist          << ' '
               << p_avg_Meas[i] << ' ' 
               << p_avg_Cond[i] << endl;
       }
   }
   os1.close();
   os2.close();

   delete [] file_1d;
   delete [] file_2d;
}

