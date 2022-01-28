#include "PressMeasure.h"

PressMeasure::PressMeasure(int num, Grid *grid, bool debug_) 
: num_p_meas(num), gridPtr(grid), debug(debug_)  
{
  if(debug) 
     cout<<"PressMeasure::PressMeasure()"
         <<endl;
  
  initilize();
}

void PressMeasure::initilize() {
    i_pmeas = new int    [num_p_meas];
    j_pmeas = new int    [num_p_meas];
    k_pmeas = new int    [num_p_meas];
  ijk_pmeas = new int    [num_p_meas];
     p_meas = new double [num_p_meas];
  
  p_avg_cond = new double [num_p_meas];
  p_avg_uncd = new double [num_p_meas];
  p_std_cond = new double [num_p_meas];
  p_std_uncd = new double [num_p_meas];
  
  w_pmeas = new double [num_p_meas];
  for(int i = 0; i < num_p_meas; ++i) {
      w_pmeas[i] = 1.;
      p_avg_cond[i] = 0.;
      p_avg_uncd[i] = 0.;
      p_std_cond[i] = 0.;
      p_std_uncd[i] = 0.;
  }
}

PressMeasure::~PressMeasure() {
  if(debug)         
     cout<<"PressMeasure::~PressMeasure()"<<endl;
  delete [] i_pmeas;
  delete [] j_pmeas;
  delete [] k_pmeas;
  delete [] ijk_pmeas;
  delete [] p_meas;
  
  delete [] p_avg_cond;
  delete [] p_avg_uncd;
  delete [] p_std_cond;
  delete [] p_std_uncd;
  
  delete [] w_pmeas;
}

void PressMeasure::setPdata(int m, int i, int j, int k, 
                            double p) {
    i_pmeas[m] = i;
    j_pmeas[m] = j;
    k_pmeas[m] = k;
  ijk_pmeas[m] = gridPtr->getIndex(i,j,k);
  p_meas[m]  = p;
}

void PressMeasure::setPuncd(int m, double p_avg) {
  p_avg_uncd[m] = p_avg;
}

void PressMeasure::setPcond(int m, double p_avg) {
  p_avg_cond[m] = p_avg;
}

void PressMeasure::setPuncd(int m, double p_avg, double p_std) {
  p_avg_uncd[m] = p_avg;
  p_std_uncd[m] = p_std;
}

void PressMeasure::setPcond(int m, double p_avg, double p_std) {
  p_avg_cond[m] = p_avg;
  p_std_cond[m] = p_std;
}

double PressMeasure::getAbsDiff(double &sumwall) {
                                
  abs_avg_diff = 0.;
  sumwall = 0;
  for(int i = 0; i < num_p_meas; ++i) {
      abs_avg_diff += fabs( p_avg_cond[i] - p_meas[i]);
      //    std_diff +=  p_std_cond[i];
      //cout << "p_std = " << p_std_cond[i]<<endl;
      sumwall   += w_pmeas[i];
  }
  abs_avg_diff /= num_p_meas;
  
  //    std_diff /= num_p_meas;
  //abs_std_diff  = std_diff; 
  return abs_avg_diff;
}

double PressMeasure::getAbsDiff_D(double &sumwall) {
                                
  abs_avg_diff = 0.;
  sumwall = 0;
  for(int i = 0; i < num_p_meas; ++i) {
      abs_avg_diff += (  fabs( p_avg_cond[i] - p_meas[i])
                       + p_std_cond[i]
                      );
      /*
      cout << p_avg_cond[i] << ' '
           << p_meas[i]     << ' '
           << p_std_cond[i] << endl;
      */
      sumwall   += w_pmeas[i];
  }
  abs_avg_diff /= num_p_meas;
  
  return abs_avg_diff;
}



double PressMeasure::getSqdDiff(){
  double sum     = 0.;
  //cout << "Ppred    :  Pmeasure" <<endl;  
  for(int i = 0; i < num_p_meas; ++i) {
      sum += w_pmeas[i] *
             (  (p_avg_cond[i] - p_meas[i]) * (p_avg_cond[i] - p_meas[i])
	      +  p_std_cond[i] * p_std_cond[i]
             );
      //cout << p_avg_cond[i] << ' ' << p_meas[i] << endl;
  }
  sum     = sum     / abs_avg_diff /abs_avg_diff;
  return sum;
}

double PressMeasure::getAbsDiff(double &std_diff,
                                double &sumwall) {
  abs_avg_diff = std_diff = 0.;
  
  sumwall = 0;
  for(int i = 0; i < num_p_meas; ++i) {
      abs_avg_diff += fabs( p_avg_cond[i] - p_meas[i]);
      //    std_diff +=  p_std_cond[i];
      //cout << "p_std = " << p_std_cond[i]<<endl;
      sumwall   += w_pmeas[i];
  }
  abs_avg_diff /= num_p_meas;
      std_diff /= num_p_meas;
  abs_std_diff  = std_diff; 
  return abs_avg_diff;
}

double PressMeasure::getSqdDiff(double &sum_std){
  double sum     = 0.;
  sum_std = 0.;
  if(debug) 
     cout << "Ppred    :  Pmeasure" <<endl;
  for(int i = 0; i < num_p_meas; ++i) {
      sum += w_pmeas[i] * ( p_avg_cond[i] - p_meas[i])
                        * ( p_avg_cond[i] - p_meas[i]);
      //sum_std += w_pmeas[i] * ( p_std_cond[i] )
      //                      * ( p_std_cond[i] );
      if(debug) 
         cout << p_avg_cond[i] << ' ' << p_meas[i] << endl;
  }
  sum     = sum     / abs_avg_diff /abs_avg_diff;
  //sum_std = sum_std / abs_std_diff /abs_std_diff;
  return sum;
}

// ===== Display and Output section =====
void PressMeasure::display() {
   cout<< "PressMeasure::display()" << endl;
   for(int i = 0; i < num_p_meas; ++i) {
       cout << i_pmeas[i]<<' '
            << j_pmeas[i]<<' '
            << k_pmeas[i]<<' '
            << p_meas[i] <<' '
            << p_avg_cond[i] << endl;
   }
}

void PressMeasure::output() {
   ofstream os("pressure_data2d.out", ios::out);
   os << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(5);

   for(int i = 0; i < num_p_meas; ++i) {
       os << i_pmeas[i]<<' '
          << j_pmeas[i]<<' '
          << k_pmeas[i]<<' '
          << p_meas[i] <<' '
          << p_avg_cond[i] << endl;
   }
   os.close();
}

void PressMeasure::output(char *file) {
   char * ending1 = "_1d.out";
   char * ending2 = "_2d.out";
   int length  = strlen(file);
   int length1 = strlen(ending1);
   int length2 = strlen(ending2);
   char *file_1d = new char [length + length1];   
   char *file_2d = new char [length + length2];
   strcpy(file_1d, file);
   strcpy(file_2d, file);
   strcat(file_1d,ending1);
   strcat(file_2d,ending2); 
   
   ofstream os1(file_1d, ios::out);
   ofstream os2(file_2d, ios::out);
   os1 << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);
   os2 << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);    

   for(int i = 0; i < num_p_meas; ++i) {
       os2 << gridPtr->getX( i_pmeas[i] ) << ' ' 
           << gridPtr->getY( j_pmeas[i] ) << ' '
           << p_meas[i] << ' ' 
           << p_avg_cond[i] << endl;
       if( i_pmeas[i] == j_pmeas[i] ) {
           double dist = sqrt(  gridPtr->getX( i_pmeas[i] ) 
                              * gridPtr->getX( i_pmeas[i] )
                              + gridPtr->getY( j_pmeas[i] )
                              * gridPtr->getY( j_pmeas[i] )
                             );
           os1 << dist          << ' '
               << p_meas[i]     << ' ' 
               << p_avg_cond[i] << endl;
       }
   }
   os1.close();
   os2.close();

   delete [] file_1d;
   delete [] file_2d;
}
