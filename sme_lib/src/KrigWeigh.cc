#include "KrigWeigh.h"

KrigWeigh::KrigWeigh(Grid* grid, Perm* perm, int num) 
 : gridPtr(grid), permPtr(perm), num_meas(num)
{  debug = false;
   if(debug) 
      cout << "KrigWeigh::KrigWeigh()" << endl;
   initialization();
}

KrigWeigh::~KrigWeigh() {
   if(debug) 
      cout << "KrigWeigh::~KrigWeigh()" << endl;
   delete []weigh;
}

void KrigWeigh::initialization() {
   length = gridPtr->getNumPermNode();
   weigh  = new double [ length * num_meas];
}

void KrigWeigh::setCovariance22(int* im, int* jm, int* km, double* cov22){ 
  for(int m2 = 0; m2 < num_meas; ++m2) {
      for(int m1 = 0; m1 < num_meas; ++m1) {
          cov22[ getIndex(m1, m2) ] = permPtr->getCorYY(im[m1], jm[m1], km[m1], 
                                                        im[m2], jm[m2], km[m2]);
          //cout<<im[m1] << jm[m1] << km[m1] 
	  //    <<im[m2] << jm[m2] << km[m2] 
	  //    <<' '<<cov22[ getIndex(m1, m2) ] <<endl;
      }
  }
  //exit(0);  
}

void KrigWeigh::setCovariance12(int* im, int* jm, int* km, double* cov12) {
   for(int k = 0; k < gridPtr->getPermNz(); ++k) {
       for(int j = 0; j < gridPtr->getPermNy(); j++) {
           for(int i = 0; i < gridPtr->getPermNx(); i++) {
               int ijk = gridPtr->getPermIndex(i,j,k);
               for(int m = 0; m < num_meas; m++) {
                   cov12[ getIndex(m, ijk) ] = permPtr->getCorYY(i, j, k,
                                                          im[m], jm[m], km[m]);
              }
           }
       }
   }
}

void KrigWeigh::calWeigh(int* im, int* jm, int* km) {
   
   double* cov22    = new double[num_meas * num_meas];
   double* cov22inv = new double[num_meas * num_meas];
   double* cov12    = new double[length   * num_meas];
   
   setCovariance22(im, jm, km, cov22);
   //cout<<"setCovariance22"<<endl;
   calCovInverse(cov22, cov22inv);
   //cout<<"setCovariance22inve"<<endl;
   setCovariance12(im, jm, km, cov12);
   //cout<<"setCovariance12"<<endl;
   double sum;
   for(int i = 0; i < length; ++i) {
       for(int j = 0; j < num_meas; ++j) {
           sum = 0.;
           for(int k = 0; k < num_meas; ++k)
               sum += cov12[ getIndex(k, i) ] * cov22inv[ getIndex(k, j) ];
           weigh[ getIndex(j, i) ] = sum;
       }
   }
   delete cov22;
   delete cov22inv;
   delete cov12;   
}

void KrigWeigh::calCovInverse(double* cov22, double* cov22inv) {
   double sum;

   double* cov22a = new double[num_meas * num_meas];
   double* cov22b = new double[num_meas * num_meas];
   for(int i = 0; i < num_meas * num_meas; ++i) cov22a[i] = cov22[i];
   for(int i = 0; i < num_meas * num_meas; ++i) cov22b[i] = 0.;
   
   int ierr = cholsky(num_meas, cov22a, cov22b);

   if(ierr > 0) {
      cout <<"ierr = "<<ierr<<endl;
      exit(0);
   }
   for(int i = 0; i < num_meas * num_meas; ++i) cov22a[i] = 0.;
   
   linv(num_meas, cov22b, cov22a);
   for(int j = 0; j < num_meas; ++j) {
       for(int i = 0; i < num_meas; ++i) {
           sum = 0.;
           for(int k = 0; k < num_meas; ++k) {
               sum += cov22a[ getIndex(k, i) ] * cov22a[ getIndex(k, j) ];
           }
           cov22inv[ getIndex (i, j) ] = sum;
       }
   }
   /*
   for(int j = 0; j < num_meas; ++j) {
       for(int i = 0; i < num_meas; ++i) {
           sum = 0;
           for(int k = 0; k < num_meas; ++k) {
               sum += cov22[ getIndex(i, k ) ] * cov22inv[ getIndex(k, j) ];
           }
           cov22a[ getIndex(i, j) ] = sum;
       }
   }*/
   delete [] cov22a;
   delete [] cov22b;
}

void KrigWeigh::calCovInverse(int meas_len, double* cov22, double* cov22inv) {
   double sum;

   double* cov22a = new double[ meas_len * meas_len ];
   double* cov22b = new double[ meas_len * meas_len ];
   for(int i = 0; i < meas_len * meas_len; ++i) cov22a[i] = cov22[i];
   for(int i = 0; i < meas_len * meas_len; ++i) cov22b[i] = 0.;
   
   int ierr = cholsky(meas_len, cov22a, cov22b);

   if(ierr > 0) {
      cout <<"ierr = "<<ierr<<endl;
      exit(0);
   }
   
   for(int i = 0; i < meas_len * meas_len; ++i) cov22a[i] = 0.;
   
   linv(meas_len, cov22b, cov22a);
   for(int j = 0; j < meas_len; ++j) {
       for(int i = 0; i < meas_len; ++i) {
           sum = 0.;
           for(int k = 0; k < meas_len; ++k) {
               sum += cov22a[ k + i * meas_len ] * cov22a[ k + j * meas_len ];
           }
           cov22inv[ i + j * meas_len ] = sum;
       }
   }
   /* check the inverse result
   for(int j = 0; j < meas_len; ++j) {
       for(int i = 0; i < meas_len; ++i) {
           sum = 0;
           for(int k = 0; k < meas_len; ++k) {
               sum += cov22[ i + k * meas_len ] * cov22inv[ k + j * meas_len ];
           }
           cov22a[ i + j * meas_len ] = sum;
           cout << cov22a[i + j * meas_len ] << ' ';
       }
       cout << endl;
   }
   */
   /*
   cout<<"KrigWeight" <<endl;
   exit(0);
   */
   delete [] cov22a;
   delete [] cov22b;
}

void KrigWeigh::condGaussDist(int cova_len, int meas_len, int* node, 
                              double *cova, double *std) {
   double* b12    = new double[cova_len * meas_len];
   double* b22    = new double[meas_len * meas_len];
   double* b22inv = new double[meas_len * meas_len];
   for(int i = 0; i < cova_len; ++i) {
       for(int j = 0; j < meas_len; j++) {
           b12[ j + i * meas_len] = cova[ node[j]  + i * cova_len ];
           //cout << b12[ j + i * meas_len] << endl;
       }
   }
   for(int i = 0; i < meas_len; ++i) {
       for(int j = 0; j < meas_len; j++) {
           b22[ j + i * meas_len] = cova[node[j]  + node[i] * cova_len];
       }
   }

   calCovInverse(meas_len, b22, b22inv);
   
   double sum;
   double tmp;
   for(int i = 0; i < cova_len; ++i) {
       sum = 0.;
       for(int k1 = 0; k1 < meas_len; k1++) {
           for(int k2 = 0; k2 < meas_len; k2++) {
               sum +=     b12[ k2 +  i * meas_len]
                      *b22inv[ k2 + k1 * meas_len]
                       *  b12[ k1 +  i * meas_len];
           }
       }
       tmp = cova[i + i * cova_len] - sum;
       if(tmp < 0) {
           std[i] = 0.;
       } else {
           std[i] = sqrt(tmp);
       }
   }              
   /*
   for(int i = 0; i < cova_len; ++i) {
       for(int j = 0; j < cova_len; ++j) {
           sum = 0.;
           for(int k1 = 0; k1 < meas_len; k1++) {
               for(int k2 = 0; k2 < meas_len; k2++) {
                   sum +=     b12[ k2 +  i * meas_len]
                          *b22inv[ k2 + k1 * meas_len]
                          *   b12[ k1 +  j * meas_len];
               }
           }
           cova[j + i * cova_len] -= sum;
       }
   } 
   */
   /*
   for(int i = 0; i < meas_len; ++i) {
       for(int j = 0; j < meas_len; j++) {
           cout << b22[ j + i * meas_len] << ' ';
       }
       cout << endl;
   }
   cout <<"its inverse = " << endl;
   for(int i = 0; i < meas_len; ++i) {
       for(int j = 0; j < meas_len; j++) {
           cout << b22inv[ j + i * meas_len] << ' ';
       }
       cout << endl;
   }
   */

   delete [] b12;
   delete [] b22;
   delete [] b22inv;
}


void KrigWeigh::output(){
   ofstream os("krig.out", ios::out);
   os << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(5);
   for(int k = 0; k < gridPtr->getPermNz(); ++k) {
       for(int j = 0; j < gridPtr->getPermNy(); j++) {
           for(int i = 0; i < gridPtr->getPermNx(); i++) {
               int ijk = gridPtr->getPermIndex(i,j,k);
               os << i << ' ' << j << ' ';
               for(int m = 0; m < num_meas; ++m)
                   os<< weigh[ getIndex(m, ijk) ]<<' ';
               os<<endl;
           } 
           os << endl;
       }
       os << endl;
   }          
   os.close();
}

void KrigWeigh::display(){
   for(int k = 0; k < gridPtr->getPermNz(); ++k) {
       for(int j = 0; j < gridPtr->getPermNy(); j++) {
           for(int i = 0; i < gridPtr->getPermNx(); i++) {
               int ijk = gridPtr->getPermIndex(i,j,k);
               for(int m = 0; m < num_meas; ++m)
                   cout<< weigh[ getIndex(m, ijk) ]<<' ';
               cout<<endl;
           } 
           cout << endl;
       }
   }          
}


//----------------------------------------------------------------------
//    Cholesky Decomposition
//    **********************
//
//    This subroutine calculates the lower triangular matrix T
//    which, when multiplied by its own transpose,
//    gives the symmetric matrix A. (from "Numerical Analysis of
//    Symmetric Matrices,"  H.R. Schwarz et al., p. 254)
//   
//
//      INPUT VARIABLES:
//
//      a(n,n)     Symmetric positive definite matrix to be
//                 decomposed (destroyed in the calculation of t)
//      t(n,n)     Lower triangular matrix solution
//      n          Dimension of the system you're decomposing
//      ierr       Error code: ierr=0 - no errors,
//                             ierr=1 - matrix a is not positive definite
//
//      NO EXTERNAL REFERENCES:
//----------------------------------------------------------------------
int KrigWeigh::cholsky(int& n, double* a, double* t) {
    int ierr = 0;
    for(int ip = 0; ip < n; ++ip) {
       if(a[ip + ip * n] <= 0.0) {
          cerr << "Warning : cholsky -- not positive definite"<<endl;
          ierr = 1;
          break;
       }
       t[ip + ip * n] = sqrt(a[ip + ip * n]);
       if( ip >= n - 1) return ierr;
       for(int k = ip + 1; k < n; ++k){
           t[k + ip * n] = a[ip + k * n]/t[ip + ip * n];
       }
       for(int i = ip + 1; i < n; ++i){
           for(int k = i; k < n; ++k){
               a[i + k * n] -= t[i + ip * n] * t[k + ip * n];
           }
       }
    }
    return ierr;
}

//-----------------------------------------------------------------------
//
//       Inverse of a Lower Triangular Matrix
//       ************************************
//
//       This subroutine finds the inverse of a lower triangular matrix A
//       and stores the answer in B. (from "Numerical Analysis of Symmetric
//       Matrices,"  H.R. Schwarz et al.,)
//
//       INPUT VARIABLES:
//
//       a(n,n)       Lower triangular matrix to be inverted
//       b(n,n)       the inverse
//       n            Dimension of the matrix you're inverting
//-----------------------------------------------------------------------
void KrigWeigh::linv(int& n, double* a, double* b){
    double sum;
    for(int i = 0; i < n; ++i) {
        if(i > 0) {
           for(int k = 0; k < i; k++) {
               sum = 0.;
               for(int j = k; j < i; ++j) {
                   sum += a[i + j * n] * b[j + k * n];
               }
               b[i + k * n] = - sum/a[i + i * n];
           } 
        }
        b[i + i * n] = 1./a[i + i * n];
    }     
}

