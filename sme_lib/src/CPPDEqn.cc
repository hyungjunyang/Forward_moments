/*
 * File: CPPDEqn.cc
 * ----------------------------------
 * Implementation for CPPDEqn class
 */
#include "CPPDEqn.h"


//===== ( CPPDEqn() ) =====
CPPDEqn::CPPDEqn( P0Eqn &P0ee, Solver &s, const Control &c, bool debug_) 
  : Eqn( P0ee ), P0e( &P0ee ), debug(debug_)  {

  if(debug) cout << "CPPDEqn::CPPDEqn()" <<endl;

  small_value = 0.00000001;

  p0  = P0e->getP0();
  initialize();
  readData();
  calcAInv( s );
  calcBs(c.getBType(), c.getBCond());

  if(debug) cout << "calcAInv is done" <<endl;
}

//===== (~CPPDEqn()) =====
CPPDEqn::~CPPDEqn() {
  if(debug) cout << "CPPDEqn::~CPPDEqn()" <<endl;

  delete[] CPP;
  delete[] CYY;
  delete[] p_std_Cond;
  delete[] p_std_UnCd;
  delete[] p_std_Meas;

  if(nx > 1) {
     delete[] b1;
     delete[] b2;
  }
  if(ny > 1) {
     delete[] b3;
     delete[] b4;
  }
  if(nz > 1) {
     delete[] b5;
     delete[] b6;
  }
  delete[] aInv;
}

void CPPDEqn::initialize(){
  CYY         = new double[ nNodePerm * nNodePerm ];
  calcCYY();
  CPP         = new double[ nNode * nNode ];
  p_std_Cond  = new double[ nNode ];
  p_std_UnCd  = new double[ nNode ];
  p_std_Meas  = new double[ nNode ];

  aInv        = new double[ nNode * nNode ];
  if(nx > 1) {
     b1       = new double[ nNode ];
     b2       = new double[ nNode ];
     for(int i = 0; i < nNode; i++)
         b1[i] = b2[i] = 0.;
  }
  if(ny > 1) {
     b3       = new double[ nNode ];
     b4       = new double[ nNode ];
     for(int i = 0; i < nNode; i++)
         b3[i] = b4[i] = 0.;
  }
  if(nz > 1) {
     b5       = new double[ nNode ];
     b6       = new double[ nNode ];
     for(int i = 0; i < nNode; i++)
         b5[i] = b6[i] = 0.;
  }
  for(int i = 0; i < nNode * nNode; i++) {
      CPP[i] = aInv[i] = 0.0;
  }
  for(int i = 0; i < nNode; i++) {
      p_std_Cond[i] = p_std_UnCd[i] = 0.;
  }
}

void CPPDEqn::readData() {
   if(debug) cout << "ReadData()" << endl;

   ifstream in_Pcondi;
   in_Pcondi.open("press_condi.in",ios::in);
   if(in_Pcondi.bad() ) {
      cerr << " Can not open file " <<  "press_condi.in"  << endl;
      exit(8);
   } else {
      if(debug) 
         cout << "File press_condi.in was opened successfully!" <<endl;
   }
   
   // masterpoints object
   long seed;
   int num_i_ms, num_j_ms, num_k_ms;
   Junk(in_Pcondi); in_Pcondi >> seed;
   if(debug) cout<<"seed = "<< seed<<endl;
   Junk(in_Pcondi); in_Pcondi >> num_i_ms >> num_j_ms >> num_k_ms;
   if(debug) cout<<"num_i_ms = " << num_i_ms << ' '<< num_j_ms << ' '
       << num_k_ms << endl;
   //num_mpts = num_i_ms * num_j_ms * num_k_ms;
   //addMasterPointObj(seed, num_i_ms, num_j_ms, num_k_ms);

   // pressure measurements object
   // int num_p_meas;
   int i_pmeas, j_pmeas, k_pmeas;
   double p_meas, p_pred;
   
   Junk(in_Pcondi); in_Pcondi >> num_p_meas;
   if(debug) cout<<"num_p_meas = " << num_p_meas <<endl;
   addPressMeasureObj(num_p_meas);
   for(int m = 0; m < num_p_meas; ++m) {
       Junk(in_Pcondi); 
       in_Pcondi >> i_pmeas>>j_pmeas>>k_pmeas>>p_meas;
       if(debug) cout << i_pmeas<<' '<< j_pmeas<<' '<< k_pmeas<<' '
                 << p_meas << endl;
       prsmPtr->setPdata(m, --i_pmeas, --j_pmeas, --k_pmeas, p_meas);
   }
   //exit(0);
}

void CPPDEqn::addPressMeasureObj(int num) {
   prsmPtr = new PressMeasure(num, &g, debug);
}

void CPPDEqn::calcCYY() {
  for(int k1 = 0; k1 < nzPerm; ++k1) {
      for(int j1 = 0; j1 < nyPerm; ++j1) {
          for(int i1 = 0; i1 < nxPerm; ++i1) {
              int ijk1 = g.getPermIndex(i1, j1, k1);
              for(int k2 = 0; k2 < nzPerm; ++k2) {
                  for(int j2 = 0; j2 < nyPerm; ++j2) {
                      for(int i2 = 0; i2 < nxPerm; ++i2) {
                          int ijk2 = g.getPermIndex(i2, j2, k2);
                          CYY[ijk2 + ijk1 * nNodePerm] = permPtr->getCYY(i2, j2, k2, 
                                                                         i1, j1, k1);
                      }
                  }
              }
          }
      }
  }
}

//===== ( calculate aInv ) =====
//define aInv[irow, jcol]
void CPPDEqn::calcAInv(Solver &s) { 
    for(int jcol = 0; jcol < nNode; ++jcol) {
        for(int j = 0; j < nNode; ++j) RHS[j] = 0.;
            RHS[jcol] = 1.0;
        s.solve( false, *this, 1, RHS );
        for(int irow = 0; irow < nNode; ++irow) {
            aInv[jcol + irow * nNode] = RHS[irow];
        }
    }
}

//===== (calcBs() ) =====
void CPPDEqn::calcBs(const BType* bType, const double *bCond) {
  int it = 0;
  if( nx > 1 ) dP0dXi[0] = P0e->getDP0DXi( 0, it ); 
  if( ny > 1 ) dP0dXi[1] = P0e->getDP0DXi( 1, it ); 
  if( nz > 1 ) dP0dXi[2] = P0e->getDP0DXi( 2, it ); 
        
  double perm_interface[3];
  if( nx > 1) {
      for(int k = 0; k < nz; k++) {
          for(int j = 0; j < ny; j++) {
              for(int i = 1; i < nx - 1; i++ ) {
                  int ijk    =   g.getIndex( i  , j, k );
                  int ijk_n1 =   g.getIndex( i-1, j, k );
                  permPtr->Perm_x( i, j, k, perm_interface );
                  b2[ijk] = - perm_interface[2] * dP0dXi[0][ijk   ]/g.getBdx( i );
                  b1[ijk] =   perm_interface[0] * dP0dXi[0][ijk_n1]/g.getBdx( i );   
              }
          }
      }
      // CONST_RATE boundary condition is a default setting!! 
      if( g.getGridType() == CELL_CEN ) {
          for(int k = 0; k < nz; k++) {
              for(int j = 0; j < ny; j++ ) {

                  int i      = 0 ;
                  int ijk    = g.getIndex( i, j, k );
                  permPtr->Perm_x( i, j, k, perm_interface );
                  b2[ijk] = - perm_interface[2] * dP0dXi[0][ijk   ]/g.getBdx( i );
  
                      i      = nx - 1;
                      ijk    = g.getIndex( i,     j, k );
                  int ijk_n1 = g.getIndex( i - 1, j, k );
                  permPtr->Perm_x( i, j, k, perm_interface );
                  b1[ijk] =   perm_interface[0] * dP0dXi[0][ijk_n1]/g.getBdx( i );
              } 
          }
      }
  }
  if( ny > 1) {
      for(int k = 0; k < nz; k++) {
          for(int j = 1; j < ny - 1; j++ ) {
              for(int i = 0; i < nx; i++) {    
                  int ijk_n1 =   g.getIndex( i, j-1, k );
                  int ijk    =   g.getIndex( i,   j, k );
                  permPtr->Perm_y( i, j, k, perm_interface );
                  b4[ijk] = - perm_interface[2] * dP0dXi[1][ijk   ]/g.getBdy( j );
                  b3[ijk] =   perm_interface[0] * dP0dXi[1][ijk_n1]/g.getBdy( j );
              }
          }
      }
      if( g.getGridType() == CELL_CEN ) {
          for(int k = 0; k < nz; k++) {
              for(int i = 0; i < nx; i++) {
                  int j = 0;
                  int ijk    =   g.getIndex( i,   j, k );
                  permPtr->Perm_y( i, j, k, perm_interface );
                  b4[ijk] = - perm_interface[2] * dP0dXi[1][ijk   ]/g.getBdy( j );
                  
                      j      = ny - 1;
                      ijk    =   g.getIndex( i,   j, k );                      
                  int ijk_n1 =   g.getIndex( i, j-1, k );
                  permPtr->Perm_y( i, j, k, perm_interface );
                  b3[ijk] =   perm_interface[0] * dP0dXi[1][ijk_n1]/g.getBdy( j );
              }
          }
      }
   }
   if(debug) cout << "bs() is done" <<endl;
}

void CPPDEqn::calcAmR(int &irow, int &jcol, double *AmB){
   for(int k  = 0; k < nz; ++k) {
       int k2 = g.toKPerm( k );
       for(int j = 0; j < ny; ++j) {
           int j2    = g.toJPerm( j );
           int j2_p1 = j2 + 1;
           int j2_n1 = j2 - 1;
           for(int i = 0; i < nx; ++i) {
               int i2    = g.toIPerm( i );
               int i2_p1 = i2 + 1;
               int i2_n1 = i2 - 1;
               int ijk = g.getIndex(i,j,k);
               int ijk2_1 = g.getPermIndex(i2_n1, j2, k2);
               int ijk2_2 = g.getPermIndex(i2_p1, j2, k2);
               int ijk2_3 = g.getPermIndex(i2, j2_n1, k2);
               int ijk2_4 = g.getPermIndex(i2, j2_p1, k2);
               for(int kc = 0; kc < nz; ++kc) {
                   int kc2 = g.toKPerm( kc );
                   for(int jc = 0; jc < ny; ++jc) {
                       int jc2    = g.toJPerm( jc );
                       int jc2_p1 = jc2 + 1;
                       int jc2_n1 = jc2 - 1;
                       for(int ic = 0; ic < nx; ++ic) {
                           int ic2    = g.toIPerm( ic );
                           int ic2_p1 = ic2 + 1;
                           int ic2_n1 = ic2 - 1;
                           int ijkc = g.getIndex(ic,jc,kc);
                           int ijkc2_1 = g.getPermIndex(ic2_n1, jc2, kc2);
                           int ijkc2_2 = g.getPermIndex(ic2_p1, jc2, kc2);
                           int ijkc2_3 = g.getPermIndex(ic2, jc2_n1, kc2);
                           int ijkc2_4 = g.getPermIndex(ic2, jc2_p1, kc2);
                           double tmp = aInv[ijk  + irow * nNode] 
                                      * aInv[ijkc + jcol * nNode];
                           AmB[ijkc2_1 + ijk2_1 * nNodePerm] += b1[ijkc] * b1[ijk] * tmp;
                           AmB[ijkc2_2 + ijk2_2 * nNodePerm] += b2[ijkc] * b2[ijk] * tmp;
                           AmB[ijkc2_3 + ijk2_3 * nNodePerm] += b3[ijkc] * b3[ijk] * tmp;
                           AmB[ijkc2_4 + ijk2_4 * nNodePerm] += b4[ijkc] * b4[ijk] * tmp;
                           
                           AmB[ijkc2_1 + ijk2_2 * nNodePerm] += b1[ijkc] * b2[ijk] * tmp;
                           AmB[ijkc2_2 + ijk2_1 * nNodePerm] += b2[ijkc] * b1[ijk] * tmp;
                           AmB[ijkc2_1 + ijk2_3 * nNodePerm] += b1[ijkc] * b3[ijk] * tmp;
                           AmB[ijkc2_3 + ijk2_1 * nNodePerm] += b3[ijkc] * b1[ijk] * tmp;
                           
                           AmB[ijkc2_1 + ijk2_4 * nNodePerm] += b1[ijkc] * b4[ijk] * tmp;
                           AmB[ijkc2_4 + ijk2_1 * nNodePerm] += b4[ijkc] * b1[ijk] * tmp;
                           AmB[ijkc2_2 + ijk2_3 * nNodePerm] += b2[ijkc] * b3[ijk] * tmp;
                           AmB[ijkc2_3 + ijk2_2 * nNodePerm] += b3[ijkc] * b2[ijk] * tmp;
                           
                           AmB[ijkc2_2 + ijk2_4 * nNodePerm] += b2[ijkc] * b4[ijk] * tmp;
                           AmB[ijkc2_4 + ijk2_2 * nNodePerm] += b4[ijkc] * b2[ijk] * tmp;
                           AmB[ijkc2_3 + ijk2_4 * nNodePerm] += b3[ijkc] * b4[ijk] * tmp;
                           AmB[ijkc2_4 + ijk2_3 * nNodePerm] += b4[ijkc] * b3[ijk] * tmp;
                       }
                   }
               }
           }
       }
   }
}

//===== ( variCondition(const Control &c) ) =====
//need CYY array !!
void CPPDEqn::variCond() {
   int icond = prsmPtr->getIpmeas(0);
   int jcond = prsmPtr->getJpmeas(0); 
   int kcond = prsmPtr->getKpmeas(0);
   int ijk_cond = g.getIndex(icond, jcond, kcond);
   if(debug) cout << icond << ' '
                  << jcond << ' '
                  << kcond << ' '
                  << ijk_cond << endl;

   double *AmR = new double[nNodePerm * nNodePerm];
   for(int i = 0; i < nNodePerm * nNodePerm; ++i) AmR[i] = 0;
   calcAmR(ijk_cond, ijk_cond, AmR);
   
   //1) counting 
   int  num_x = 0;
   for(int k1 = 0; k1 < nNodePerm; ++k1) {
       if(fabs(AmR[k1 + k1 * nNodePerm]) > small_value ) {
          num_x += 1;
       }
       for(int k2 = k1 + 1; k2 < nNodePerm; ++k2) {
           if(fabs(AmR[k2 + k1 * nNodePerm]) > small_value ) {
              num_x += 1;             
           }
       }
   }
   //2) bookkeeping
   int    *k1_x  = new int   [num_x];
   int    *k2_x  = new int   [num_x];
   double *AmM1  = new double[num_x];
   double *AmM3  = new double[num_x];
   double *AmO   = new double[num_x];
   double *xSolu = new double[num_x];
   num_x = 0;
   for(int k1 = 0; k1 < nNodePerm; ++k1) {
       if(fabs(AmR[k1 + k1 * nNodePerm]) > small_value) {
          k1_x[num_x] = k1;
          k2_x[num_x] = k1;
          AmM3[num_x] = AmR[k1 + k1 * nNodePerm];
          AmM1[num_x] = CYY[k1 + k1 * nNodePerm];
          num_x += 1;              
       }
       for(int k2 = k1 + 1; k2 < nNodePerm; ++k2) {
           if(fabs(AmR[k2 + k1 * nNodePerm]) > small_value) {
              k1_x[num_x] = k1;
              k2_x[num_x] = k2;
              AmM3[num_x] = AmR[k2 + k1 * nNodePerm] * 2.;
              AmM1[num_x] = CYY[k2 + k1 * nNodePerm];
              num_x += 1;                       
           }
       }
   }
   /*
   double sum = 0;
   for(int i = 0; i < num_x; ++i) {
       sum += AmM1[i] * AmM3[i];
   }
   cout << sum << ' ' << sqrt(sum) << endl;
   exit(0);
   */
   // 3) calculate AmO 
   //int ijk_amo = g.getIndex(icond - 1, jcond - 1, kcond);
   int ijk_amo = ijk_max;
   //cout << g.getIndex(icond - 1, jcond - 1, kcond) << ' '
   //        << ijk_max << endl;
   //exit(0);
       
   for(int i = 0; i < nNodePerm * nNodePerm; ++i) AmR[i] = 0;
   calcAmR(ijk_amo, ijk_amo, AmR);
   for(int i = 0; i < num_x; ++i) {
       AmO[i] = AmR[k2_x[i] + k1_x[i] * nNodePerm];
   }
   delete [] AmR;
   if(debug) cout << "start optimization! " << num_x << endl;
   calcPermStdBySimplx(num_x, AmM1, AmM3, AmO, xSolu);
   if(debug) cout << "end optimization!" << endl;
   for(int i = 0; i < num_x; ++i) {
       CYY[k2_x[i] + k1_x[i] * nNodePerm] = xSolu[i];
       CYY[k1_x[i] + k2_x[i] * nNodePerm] = xSolu[i];
       if(xSolu[i] > AmM1[i])
          cout << i << ' ' 
               << "xSolu[i] = " << xSolu[i] << " > "
               << "AmM1[i] = " << AmM1[i] << endl; 
   }
   
   delete[] k1_x;
   delete[] k2_x;
   delete[] AmM1;
   delete[] AmM3;
   delete[] AmO;
   delete[] xSolu;

   // after conditioning
   CPPSolution();
   
   YVarInterpolation(CYY);
   printCYY(CYY);
}

//===== ( CPPSolution ) =====
// need CYY array !! 
void CPPDEqn::CPPSolution() {
   for(int jcol = 0; jcol < nNode; ++jcol) {
       for(int irow = jcol; irow < nNode; ++irow) {
           double sum = 0.;
           for(int k  = 0; k < nz; ++k) {
               int k2 = g.toKPerm( k );
               for(int j = 0; j < ny; ++j) {
                   int j2    = g.toJPerm( j );
                   int j2_p1 = j2 + 1;
                   int j2_n1 = j2 - 1;
                   for(int i = 0; i < nx; ++i) {
                       int i2    = g.toIPerm( i );
                       int i2_p1 = i2 + 1;
                       int i2_n1 = i2 - 1;
                       int ijk = g.getIndex(i,j,k);
                       int ijk2_1 = g.getPermIndex(i2_n1, j2, k2);
                       int ijk2_2 = g.getPermIndex(i2_p1, j2, k2);
                       int ijk2_3 = g.getPermIndex(i2, j2_n1, k2);
                       int ijk2_4 = g.getPermIndex(i2, j2_p1, k2);
                       for(int kc = 0; kc < nz; ++kc) {
                           int kc2 = g.toKPerm( kc );
                           for(int jc = 0; jc < ny; ++jc) {
                               int jc2    = g.toJPerm( jc );
                               int jc2_p1 = jc2 + 1;
                               int jc2_n1 = jc2 - 1;
                               for(int ic = 0; ic < nx; ++ic) {
                                   int ic2    = g.toIPerm( ic );
                                   int ic2_p1 = ic2 + 1;
                                   int ic2_n1 = ic2 - 1;
                                   int ijkc = g.getIndex(ic,jc,kc);
                                   int ijkc2_1 = g.getPermIndex(ic2_n1, jc2, kc2);
                                   int ijkc2_2 = g.getPermIndex(ic2_p1, jc2, kc2);
                                   int ijkc2_3 = g.getPermIndex(ic2, jc2_n1, kc2);
                                   int ijkc2_4 = g.getPermIndex(ic2, jc2_p1, kc2);
                                   sum += aInv[ijk  + irow * nNode]
                                        * aInv[ijkc + jcol * nNode]
                                        * ( b1[ijkc] * b1[ijk] * CYY[ijkc2_1 + ijk2_1 * nNodePerm]
                                          + b2[ijkc] * b2[ijk] * CYY[ijkc2_2 + ijk2_2 * nNodePerm]
                                          + b3[ijkc] * b3[ijk] * CYY[ijkc2_3 + ijk2_3 * nNodePerm]
                                          + b4[ijkc] * b4[ijk] * CYY[ijkc2_4 + ijk2_4 * nNodePerm]
                                          + b1[ijkc] * b2[ijk] * CYY[ijkc2_1 + ijk2_2 * nNodePerm]
                                          + b2[ijkc] * b1[ijk] * CYY[ijkc2_2 + ijk2_1 * nNodePerm]
                                          + b1[ijkc] * b3[ijk] * CYY[ijkc2_1 + ijk2_3 * nNodePerm]
                                          + b3[ijkc] * b1[ijk] * CYY[ijkc2_3 + ijk2_1 * nNodePerm]
                                          + b1[ijkc] * b4[ijk] * CYY[ijkc2_1 + ijk2_4 * nNodePerm]
                                          + b4[ijkc] * b1[ijk] * CYY[ijkc2_4 + ijk2_1 * nNodePerm]
                                          + b2[ijkc] * b3[ijk] * CYY[ijkc2_2 + ijk2_3 * nNodePerm]
                                          + b3[ijkc] * b2[ijk] * CYY[ijkc2_3 + ijk2_2 * nNodePerm]
                                          + b2[ijkc] * b4[ijk] * CYY[ijkc2_2 + ijk2_4 * nNodePerm]
                                          + b4[ijkc] * b2[ijk] * CYY[ijkc2_4 + ijk2_2 * nNodePerm]
                                          + b3[ijkc] * b4[ijk] * CYY[ijkc2_3 + ijk2_4 * nNodePerm]
                                          + b4[ijkc] * b3[ijk] * CYY[ijkc2_4 + ijk2_3 * nNodePerm]
                                          );
                                }
                           }
                       }
                   }
               }
           }
           CPP[jcol + irow * nNode] = sum;
       }
   }
   for(int jcol = 0; jcol < nNode; ++jcol) {
       for(int irow = 0; irow < jcol; ++irow) {
           CPP[jcol + irow * nNode] = CPP[irow  + jcol * nNode];
       }
   }   
}

void CPPDEqn::calcCondPStd() {
   for(int i = 0; i < nNode; ++i) {
       if( CPP[i + i * nNode] < 1.e-10 ) {
            p_std_Cond[i] = 0.;
       } else {
            p_std_Cond[i] = sqrt( CPP[i + i * nNode] );
       }
   }
   if(debug) display(p_std_Cond);
   printCPP("Cond_CPP.out");
}

void CPPDEqn::calcUnCdPStd() {
   for(int i = 0; i < nNode; ++i) {
       if( CPP[i + i * nNode] < 1.e-10 ) {
            p_std_UnCd[i] = 0.;
       } else {
            p_std_UnCd[i] = sqrt( CPP[i + i * nNode] );
       }
   }
   if(debug) display(p_std_UnCd);
   printCPP("UnCd_CPP.out");
}

void CPPDEqn::calcPermStdBySimplx(int& num_x, double *AmM1, double *AmM3, double *AmO, double *xSolu) {
  int icase;
  int n  = num_x;
  int np = n + 1;   
  int m1 = n;
  int m2 = 0; 
  int m3 = 1;
  int m  = m1 + m2 + m3; 
  int mp = m + 2;
  int nm1m2 = n + m1 + m2;
  int m1p1 = m1 + 1;
  cout << "the size of am = " << mp * np << endl;
  double * am    = new double[mp * np];
  int    * izrov = new int [n];
  int    * iposv = new int [m];

  // Initialize to 0.
  for(int i = 0; i < mp * np; ++i) am[i] = 0.0;
  
  // The first row in the tableau, or row 0
  am[ 0 + 0 * mp] =  0.0;
  for(int j = 1; j < np; ++j) {
      //am[ 0 + j * mp] =  1.0;
      am[ 0 + j * mp] =  AmO[j - 1];
  }
 
  // Rows from 1 to m1 in the tableau, the part with <=  
  for(int i = 1; i <= m1; ++i) {
    //am[ i + 0 * mp] =  1.0;
      am[ i + 0 * mp] =  AmM1[ i - 1 ];
      am[ i + i * mp] = -1.0;
    //am[ i + i * mp] = -AmM1[ i - 1 ];  
  }
  for(int i = 0; i < m3; ++i) {
          am[ m1p1 + i + 0 * mp] = 0.0;
      for(int j = 1; j < np; ++j) {
          am[ m1p1 + i + j * mp] = AmM3[j - 1];
      }
  }
  int nmax = n;
  int mmax = m;
  int    * l1    = new int [nmax];
  int    * l2    = new int [mmax];
  int    * l3    = new int [mmax];
  simplx(am, &m, &n, &mp, &np,
         &m1, &m2, &m3, &icase, izrov, iposv,
         &nmax, l1, &mmax, l2, l3);
  //cout << icase << endl;
  //exit(0);
  if(        icase ==  1) {
     cout << "Unbounded objective function " << endl;
  } else if (icase == -1) {
     cout << "No solution satisfy constraints given " << endl;
  } else {
     for(int i = 0; i < m; ++i) {
         if(iposv[i] <= n) {
            xSolu[iposv[i] - 1] = am[ i  + 1];
        }
     }
  }
  delete [] izrov;
  delete [] iposv;
  delete [] am;
  delete [] l1;
  delete [] l2;
  delete [] l3;
}

//===== ( Solve ) =====
void CPPDEqn::solve() {
   
   double b1b1, b2b2, b3b3, b4b4;
   double b1b2, b2b1, b1b3, b3b1;
   double b1b4, b4b1, b2b3, b3b2;
   double b2b4, b4b2, b3b4, b4b3;

   for(int jcol = 0; jcol < nNode; ++jcol) {
       for(int irow = 0; irow < nNode; ++irow) {
           double sum = 0.;
           for(int k  = 0; k < nz; ++k) {
               int k2 = g.toKPerm( k );
               for(int j = 0; j < ny; ++j) {
                   int j2    = g.toJPerm( j );
                   int j2_p1 = j2 + 1;
                   int j2_n1 = j2 - 1;
                   for(int i = 0; i < nx; ++i) {
                       int i2    = g.toIPerm( i );
                       int i2_p1 = i2 + 1;
                       int i2_n1 = i2 - 1;
                       int ijk = g.getIndex(i,j,k);
                       for(int kc = 0; kc < nz; ++kc) {
                           int kc2 = g.toKPerm( kc );
                           for(int jc = 0; jc < ny; ++jc) {
                               int jc2    = g.toJPerm( jc );
                               int jc2_p1 = jc2 + 1;
                               int jc2_n1 = jc2 - 1;
                               for(int ic = 0; ic < nx; ++ic) {
                                   int ic2    = g.toIPerm( ic );
                                   int ic2_p1 = ic2 + 1;
                                   int ic2_n1 = ic2 - 1;
                                   int ijkc = g.getIndex(ic,jc,kc);
                                   b1b1 = permPtr->getCYY(i2_n1,  j2,  k2,
                                                         ic2_n1, jc2, kc2);
                                   b2b2 = permPtr->getCYY(i2_p1,  j2,  k2,
                                                         ic2_p1, jc2, kc2);
                                   b3b3 = permPtr->getCYY(i2, j2_n1,  k2,
                                                         ic2,jc2_n1, kc2);
                                   b4b4 = permPtr->getCYY(i2, j2_p1,  k2,
                                                         ic2,jc2_p1, kc2);
                                   b1b2 = permPtr->getCYY(i2_n1,  j2,  k2,
                                                         ic2_p1, jc2, kc2);
                                   b2b1 = permPtr->getCYY(i2_p1,  j2,  k2,
                                                          ic2_n1, jc2, kc2);
                                   b1b3 = permPtr->getCYY(i2_n1,  j2,  k2,
                                                         ic2, jc2_n1, kc2);
                                   b3b1 = permPtr->getCYY(i2,    j2_n1, k2,
                                                         ic2_n1,jc2,   kc2);
                                   b1b4 = permPtr->getCYY(i2_n1,  j2,  k2,
                                                          ic2,    jc2_p1,kc2);
                                   b4b1 = permPtr->getCYY(i2,     j2_p1, k2,
                                                         ic2_n1, jc2, kc2);
                                   b2b3 = permPtr->getCYY(i2_p1,  j2,     k2,
                                                          ic2,   jc2_n1, kc2);
                                   b3b2 = permPtr->getCYY(i2,    j2_n1, k2,
                                                         ic2_p1,jc2,   kc2);
                                   b2b4 = permPtr->getCYY(i2_p1,  j2,  k2,
                                                         ic2,    jc2_p1,kc2);
                                   b4b2 = permPtr->getCYY(i2,    j2_p1,  k2,
                                                         ic2_p1,jc2,   kc2);
                                   b3b4 = permPtr->getCYY(i2,  j2_n1,  k2,
                                                         ic2, jc2_p1, kc2);
                                   b4b3 = permPtr->getCYY(i2,  j2_p1,  k2,
                                                         ic2, jc2_n1, kc2);
                                   sum += aInv[ijk  + irow * nNode]
                                        * aInv[ijkc + jcol * nNode]
                                        * ( b1[ijk] * b1[ijkc] * b1b1
                                          + b2[ijk] * b2[ijkc] * b2b2
                                          + b3[ijk] * b3[ijkc] * b3b3
                                          + b4[ijk] * b4[ijkc] * b4b4
                                          + b1[ijk] * b2[ijkc] * b1b2
                                          + b2[ijk] * b1[ijkc] * b2b1
                                          + b1[ijk] * b3[ijkc] * b1b3
                                          + b3[ijk] * b1[ijkc] * b3b1
                                          + b1[ijk] * b4[ijkc] * b1b4
                                          + b4[ijk] * b1[ijkc] * b4b1
                                          + b2[ijk] * b3[ijkc] * b2b3
                                          + b3[ijk] * b2[ijkc] * b3b2
                                          + b2[ijk] * b4[ijkc] * b2b4
                                          + b4[ijk] * b2[ijkc] * b4b2
                                          + b3[ijk] * b4[ijkc] * b3b4
                                          + b4[ijk] * b3[ijkc] * b4b3
                                          );
                                }
                           }
                       }
                   }
               }
           }
           CPP[jcol + irow * nNode] = sum;
       }
   }
   printCPP("UnCd_CPP.out");
}

void CPPDEqn::calcCPPMeasbySolveEqn(int nWells, const Well *w, 
                                    const Control &c, Solver &s, CYPEqn *CYPee, 
                                    CPPEqn *CPPee) {
   // condition part
   P0e->solve(nWells, w, c, s, 0, num_p_meas, 
              prsmPtr->getIpmeas(),
              prsmPtr->getJpmeas(),
              prsmPtr->getKpmeas(),
              prsmPtr->getPmeas()
             );
   CYPee->solve(nWells, w, c, s, 0);
   CPPee->solve(nWells, w, c, s, 0);

   for(int i = 0; i < nNode; ++i) {
       if( CPPee->getCPP(i,i) < 1.e-10 ) {
            p_std_Meas[i] = 0.;
       } else {
            p_std_Meas[i] = sqrt( CPPee->getCPP(i,i));
       }
   }
   outputObjCondCPP();
   
   double max = 0;
   for(int k = 0; k < g.getNz(); k++) {
       for(int j = 0; j < prsmPtr->getJpmeas(0); j++) {
           int i = j;
           int ijk = g.getIndex(i, j, k);
           if(p_std_Meas[ ijk ] > max) {
              ijk_max = ijk;
              max = p_std_Meas[ ijk ];
           }
           //cout << ijk << ' ' << p_std_Meas[ ijk ] << endl;
       }
   } 
}

/*
void CPPDEqn::condCPPbyGaussDist() {
   
   krigPtr->condGaussDist(gridPtr->getnNode(), num_p_meas, 
                          prsmPtr->getIJKpmeas(), CPPe->getCPP(),
                          P_std_uncd_meas
                          );
   
   outputObjCondCPP();
}
*/

void CPPDEqn::YVarInterpolation(double *CYY) {
   for(int k = 0; k < g.getNz(); k++) {
       int k2 = g.toKPerm( k );
       for(int j = 0; j < g.getNy(); j++) {
           int j2 = g.toJPerm( j );
	   int j2_p1 = j2 + 1;
           int j2_n1 = j2 - 1;
           for(int i = 0; i < g.getNx(); i++) {
               int i2    = g.toIPerm( i );
               int i2_p1 = i2 + 1;
               int i2_n1 = i2 - 1;
	       int ijk2     = g.getPermIndex(i2   , j2   , k2);
	       int ijk2_ip1 = g.getPermIndex(i2_p1, j2   , k2);
               int ijk2_in1 = g.getPermIndex(i2_n1, j2   , k2);
 	       int ijk2_jp1 = g.getPermIndex(i2   , j2_p1, k2);
               int ijk2_jn1 = g.getPermIndex(i2   , j2_n1, k2);
	       CYY[ijk2 + ijk2 * nNodePerm] = (  CYY[ijk2_ip1 + ijk2_ip1 * nNodePerm]
                                               + CYY[ijk2_in1 + ijk2_in1 * nNodePerm]
                                               + CYY[ijk2_jp1 + ijk2_jp1 * nNodePerm]
                                               + CYY[ijk2_jn1 + ijk2_jn1 * nNodePerm]
			                      )/4.;
           }
       } 
   }          

}

//===== output() =====
void CPPDEqn::printCYY(double *CYY) {
   ofstream os2("CondiPerm_2d.out", ios::out);
   os2 << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(5);
   for(int k = 0; k < g.getPermNz(); k++) {
       for(int j = 0; j < g.getPermNy(); j++) {
           for(int i = 0; i < g.getPermNx(); i++) {
               int ijk = g.getPermIndex(i, j, k);
               os2 << i << ' ' 
                   << j << ' '
                   << CYY[ijk + ijk * nNodePerm] <<' '
                   << endl;
           }
           os2<<endl;
       } 
   }          
   os2<<endl;
   os2.close();
}

void CPPDEqn::printCPP(char *file){
   ofstream os(file, ios::out);
   for(int j = 0; j < nNode; j++) 
       os << CPP[j + j * nNode] << endl;
   os<<endl;
   for(int j = 0; j < nNode; j++) {
       for(int i = 0; i < nNode; i++) {
           os << CPP[i + j * nNode] << endl;
       }
   }
   os.close();
}

void CPPDEqn::output(){
   ofstream os("CondiPressure_2d.out", ios::out);
   os << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(5);
   for(int k = 0; k < g.getNz(); k++) {
       for(int j = 0; j < g.getNy(); j++) {
           for(int i = 0; i < g.getNx(); i++) {
               int ijk = g.getIndex(i, j, k);
               os << g.getX(i) << ' ' 
                  << g.getY(j) << ' '
                  << p_std_UnCd[ijk]
                  << p_std_Cond[ijk]  
                  << endl;
           }
           os<<endl;
       } 
   }          
   os<<endl;
   os.close();
   
   ofstream os1("CondiPressure_1d.out", ios::out);
   os1 << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(5);
   cout << "print CondiPressure_1d.out is invoked!" << endl;
   for(int k = 0; k < g.getNz(); k++) {
       for(int j = 0; j < g.getNy(); j++) {
           int i = j;
           int ijk = g.getIndex(i, j, k);
           double dist = sqrt(  g.getX(i) * g.getX(i) 
                              + g.getY(i) * g.getY(i)
                             );
           os1 << dist   << ' '
               << p_std_UnCd[ijk]<<' '
               << p_std_Cond[ijk] 
               << endl;
       } 
   }          
   os1<<endl;
   os1.close();
}

void CPPDEqn::outputObjCondCPP() {
   ofstream os("ObjCondStdP_1d.out", ios::out);
   os << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(5);
   for(int k = 0; k < g.getNz(); k++) {
       for(int j = 0; j < g.getNy(); j++) {
           int i = j;
           int ijk = g.getIndex(i, j, k);
           double dist = sqrt(  g.getX(i) * g.getX(i) 
                              + g.getY(i) * g.getY(i) );
           os << dist    << ' '
              << p0[ijk] << ' '
              << p_std_Meas[ ijk ] << endl;
       } 
   }          
   os<<endl;
   os.close();
}

void CPPDEqn::display(double *p_std) {
   cout << setiosflags(ios::fixed | ios::showpoint)
        << setprecision(5);

   for(int irow = ny - 1; irow >= 0; --irow) {
       for(int jcol = 0; jcol < nx; ++jcol) {
           int ijk = jcol + irow * nx;
           cout << p_std[ijk] << ' ';
       }
       cout << endl;
   }
}

/* --- Private Methods ----------------------------- */
// --- C++ requirement to call Fortran subroutine ----
extern "C"
{ 
  void simplx_(double *amatrix, int *m, int *n, int *mp, int *np,
               int *m1, int *m2, int *m3, int *icase, int *izrov, int *iposv,
               int *nmax, int *l1, int* mmax, int *l2, int *l3);
};

void CPPDEqn::simplx(double *amatrix, int *m, int *n, int *mp, int *np, int *m1, 
                     int *m2, int *m3, int *icase, int *izrov, int *iposv,
                     int *nmax, int *l1, int* mmax, int *l2, int *l3) {
  simplx_(amatrix, m, n, mp, np, m1, m2, m3, icase, izrov, iposv,
          nmax, l1, mmax, l2, l3);
  if(debug) cout << "simplx is done" << endl;
}

