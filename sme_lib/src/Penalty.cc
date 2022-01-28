#include "Penalty.h"

Penalty::Penalty(Perm* perm, MasterPoint* mpts, PressMeasure *prsm, 
                 P0_Sensitivity *sens) 
: permPtr(perm), mptsPtr(mpts), prsmPtr(prsm), sensPtr(sens)
{
   cout<<"Penalty::Penalty()"<<endl;
   initilize();
}

Penalty::~Penalty() {
   cout<<"Penalty::~Penalty()"<<endl;   
   delete [] cpq;
   delete [] ppq;
   delete [] xpq;
   delete [] one_neg;
   delete [] one_pos;
   delete [] bpq_min;
   delete [] bpq_max;

   delete [] Y_avg;
   delete [] Y_var;
   delete [] Y_avg_cond;
   delete [] Y_var_cond;
   delete [] index_active;
   delete [] grad;
   delete [] u;
   delete [] uprim;
   delete [] hh_t1;
      	   
   delete [] p_pred_meas;
   delete [] p_pred_meas_wpms;
}

void Penalty::initilize() {
   num_mpts = mptsPtr->getLength();        
   cpq     = new double [ num_mpts * num_mpts ];
   ppq     = new double [ num_mpts ];
   xpq     = new double [ num_mpts ];
   one_neg = new double [ num_mpts ];
   one_pos = new double [ num_mpts ];
   bpq_min = new double [ num_mpts ]; 
   bpq_max = new double [ num_mpts ]; 
    
   Y_avg        = new double [ num_mpts ];
   Y_var        = new double [ num_mpts ];
   Y_avg_cond   = new double [ num_mpts ];
   Y_var_cond   = new double [ num_mpts ];
   grad         = new double [ num_mpts ];
   index_active = new int    [ num_mpts ];
   u            = new double [ num_mpts ];
   uprim        = new double [ num_mpts ];
   hh_t1        = new double [ num_mpts ];
   for(int i = 0; i < num_mpts; ++i)
       grad[i] = u[i] = uprim[i] = hh_t1[i] = 0.;

   num_p_meas       = prsmPtr->getLength();
   p_pred_meas      = new double [ num_p_meas ];
   p_pred_meas_wpms = new double [ num_p_meas ];
}

void Penalty::setY0data() {
   for(int mst = 0; mst < num_mpts; ++mst ) {
       Y_avg[mst] = permPtr->getYAvg( mptsPtr->getIm(mst), 
                                      mptsPtr->getJm(mst), 
                                      mptsPtr->getKm(mst) );
       Y_var[mst] = permPtr->getYVar( mptsPtr->getIm(mst), 
                                      mptsPtr->getJm(mst), 
                                      mptsPtr->getKm(mst) );
   }
}

void Penalty::setYConddata() {
   for(int mst = 0; mst < num_mpts; ++mst ) {
       Y_avg_cond[mst] = permPtr->getYAvg( mptsPtr->getIm(mst), 
                                           mptsPtr->getJm(mst), 
                                           mptsPtr->getKm(mst) );
       Y_var_cond[mst] = permPtr->getYVar( mptsPtr->getIm(mst), 
                                           mptsPtr->getJm(mst), 
                                           mptsPtr->getKm(mst) );
   }
}


void Penalty::setP0data() {
   for(int i = 0; i < num_p_meas; ++i) {
      p_pred_meas[ i ]      = prsmPtr->getPmeas(i) - prsmPtr->getPpred(i);
      p_pred_meas_wpms[ i ] = p_pred_meas[ i ] * prsmPtr->getWpmeas(i);
   }
}

void Penalty::initObjValue() {
  avrg_diff = 0.;
  for(int i = 0; i < num_p_meas; ++i) {
      avrg_diff += fabs( p_pred_meas[i]);
  }
  avrg_diff /= num_p_meas;
  avrg_diff2 = avrg_diff * avrg_diff;

  calcObjValue();
  fobj_init = fobj_iter;
  fobj_norm = 1.;
  output();
}

void Penalty::calcObjValue() {
  fobj_iter = 0.;
  for(int i = 0; i < num_p_meas; ++i) 
      fobj_iter += p_pred_meas_wpms[ i ]/avrg_diff 
                      * p_pred_meas[ i ]/avrg_diff;
}

void Penalty::iterObjValue() {
  calcObjValue();
  fobj_norm = fobj_iter/ fobj_init;
  output();
}

void Penalty::buildMatrix() {
  double ac;
  // calculate cpq
  for(int mst2 = 0; mst2 < num_mpts; ++mst2) {
      for(int mst1 = mst2; mst1 < num_mpts; ++mst1) {
          ac = 0.;
          for(int msm = 0; msm < num_p_meas; ++msm ) {
              ac +=   prsmPtr->getWpmeas( msm ) 
                    * sensPtr->getP0Sensitivity(msm, mst1)
                    * sensPtr->getP0Sensitivity(msm, mst2);
          }
          cpq[ mst1 + mst2 * num_mpts ] = ac / avrg_diff2;
          cpq[ mst2 + mst1 * num_mpts ] = ac / avrg_diff2;
      }
  } 
  // calculate ppq
  for(int mst = 0; mst < num_mpts; ++mst) {
      ac = 0.;
      for(int msm = 0; msm < num_p_meas; ++msm ) {
          ac +=  p_pred_meas_wpms[ msm ] 
               * sensPtr->getP0Sensitivity(msm, mst);
      }
      ppq[ mst ] = 2. * ac / avrg_diff2;
  }
  // one_neg, one_pos, bpq_min, bpq_max
  double Y_min, Y_max;
  for(int mst = 0; mst < num_mpts; ++mst) {
      one_neg[mst] = -1.;
      one_pos[mst] =  1.;
      Y_min  = Y_avg[mst] - 2.0 * sqrt( Y_var[mst] ); 
      Y_max  = Y_avg[mst] + 2.0 * sqrt( Y_var[mst] );
      if( Y_avg_cond[mst] > Y_max) Y_max = Y_avg_cond[mst];
      if( Y_avg_cond[mst] < Y_min) Y_min = Y_avg_cond[mst];
      bpq_min[mst] = Y_avg_cond[mst] - Y_min;
      bpq_max[mst] = Y_max - Y_avg_cond[mst];
      xpq[mst] = ( bpq_min[mst] - bpq_max[mst])/2./one_neg[mst];
  }
}

void Penalty::calcGrad() {
  double ac;
  for(int mst2 = 0; mst2 < num_mpts; ++mst2) {
      ac = 0.;
      for(int mst1 = 0; mst1 < num_mpts; ++mst1) {
          ac +=   cpq[ mst1 + mst2 * num_mpts ] * xpq[mst1];
      }
      grad[mst2] = ac * 2. + ppq[ mst2 ];
  } 
}

void Penalty::calcNumActRst() {
   num_active = 0;
   double df_min, df_max;
   for(int mst = 0; mst < num_mpts; ++mst) {
       df_min = -  one_neg[mst] * xpq[mst] + bpq_min[mst];
       df_max = -  one_pos[mst] * xpq[mst] + bpq_max[mst];
       if        ( df_min < 0.) {
           index_active[num_active] = mst;
           ++num_active;
           xpq[mst] = bpq_min[mst]/one_neg[mst];
       } else if ( df_max < 0.) {
	   index_active[num_active] = mst;
           ++num_active;
           xpq[mst] = bpq_max[mst]/one_pos[mst];
       }
   }
}

void Penalty::projectGrad() {
   for(int i = 0; i < num_active; ++i) {
       uprim[i] = 0.0;
       hh_t1[i] = 1.0;
   }
   for(int i = 0; i < num_active; ++i) {
       uprim[i] = 0.0;
       hh_t1[i] = 1.0;
   }
}

void Penalty::minimize() {
    calcGrad();
    calcNumActRst();
    for(int mst = 0; mst < num_mpts; ++mst) 
	u[ mst ] = - grad[ mst ];
    if( num_active != 0 ) projectGrad();
}

// ===== output() =====
void Penalty::output() {
  cout << "The normalized Penalty = " << fobj_norm << endl;
  ofstream os("penalty.out", ios::out);
  os << setiosflags(ios::fixed | ios::showpoint)
     << setprecision(5);
  for(int mst2 = 0; mst2 < num_mpts; ++mst2) {
      for(int mst1 = 0; mst1 < num_mpts; ++mst1) {
          os << cpq[ mst1 + mst2 * num_mpts ] << ' ';
      }
      os << ppq[mst2] <<' '
         << bpq_min[mst2]<<' '
         << xpq[mst2] <<' '
         << bpq_max[mst2]<<' '
         <<endl;
  } 
  os.close();
}

