#include "MasterPoint.h"

using namespace std;

MasterPoint::MasterPoint(long iseed, int num_i, int num_j, int num_k)
 : seed(iseed), num_i_ms(num_i), num_j_ms(num_j), num_k_ms(num_k){
   debug = false;
   if(debug) cout << "MasterPoint::MasterPoint()" << endl;
   num_ms = num_i_ms * num_j_ms * num_k_ms;
   initialization();
}

MasterPoint::~MasterPoint(){
   if(debug) cout << "MasterPoint::~MasterPoint()" << endl;
   delete[] i_ms;
   delete[] j_ms;
   delete[] k_ms;
   delete[] ijk_ms;
   delete[] x_ms;
   delete[] y_ms;
   delete[] z_ms;
}

void MasterPoint::initialization(){
     i_ms = new    int[num_ms];
     j_ms = new    int[num_ms];
     k_ms = new    int[num_ms];
   ijk_ms = new    int[num_ms];   
     x_ms = new double[num_ms];
     y_ms = new double[num_ms];
     z_ms = new double[num_ms];
}

void MasterPoint::calMasterPts(int nx, int ny, int nz, 
		               double* x1, double* x2, double* x3) {
   int ijk;
   int imst , jmst , kmst;
   int JZone, IZone, KZone;
   int j1, j2, i1, i2, k1, k2;
   double rando1, rando2, rando3;

   for(int k = 0; k < num_k_ms; ++k) {
       KZone =  nz / (num_k_ms );
       k1 = KZone * k;
       k2 = KZone * (k + 1) - 1;
       if(k == num_k_ms - 1) k2 = nz - 1;
       KZone = k2 - k1 + 1;
       //cout << k1 <<' ' <<k2 <<' '<< KZone << endl;

       for(int j = 0; j < num_j_ms; ++j) {
           JZone =  ny / (num_j_ms );
           j1 = JZone * j;
           j2 = JZone * (j + 1) - 1;
           if(j == num_j_ms - 1) j2 = ny - 1;
           JZone = j2 - j1 + 1;
       
           for(int i = 0; i < num_i_ms; ++i) {
               IZone =  nx / (num_i_ms );
               i1 = IZone * i;
               i2 = IZone * (i+1) - 1;
               if(i == num_i_ms - 1) i2 = nx - 1;
               IZone = i2 - i1 + 1;
	   
               ijk = getIndex(i,j, k);
	   
               bool repeat = false;
               do {
		  do { 
                     rando1 = ran3(&seed);
                     kmst   = int( rando1 * KZone)  + k1; 
		     //cout <<"kmst = "<< kmst << endl;
                  } while (kmst < k1 || kmst > k2);
	          //cout << "kmst out = "<< kmst << endl;
		  //exit(0);
                  do { 
                     rando2 = ran3(&seed);
                     jmst   = int( rando2 * JZone)  + j1; 
                  } while (jmst < j1 || jmst > j2);
		  
                  do {
                     rando3 = ran3(&seed);
                     imst   = int(rando3 * IZone)  + i1;
                  } while (imst < i1 || imst > i2);

                  for(int k = 0; k < ijk; ++k) {
                      if(imst == i_ms[k] && jmst == j_ms[k] ) {
                         repeat = true;
                      }
                  }
               } while( repeat);

               i_ms[ijk] =    imst;
	       j_ms[ijk] =    jmst;
	       k_ms[ijk] =    kmst;

	       ijk_ms[ijk] = imst + nx * jmst + ny * nx * kmst;
	       
               x_ms[ijk] = x1[imst];
               y_ms[ijk] = x2[jmst];
	       z_ms[ijk] = x3[kmst];
	   }
       }
   }
   
   //cout <<"done in calculating master points"<< endl;
}

void MasterPoint::setMstPts(double* x1, double* x2, double* x3) {

   /*		
   i_ms[0] = 80;
   j_ms[0] = 79;
   i_ms[1] = 78;
   j_ms[1] = 79;
   i_ms[2] = 79;
   j_ms[2] = 80;
   i_ms[3] = 79;
   j_ms[3] = 78;
   */
   
   i_ms[0] = 40;
   j_ms[0] = 39;
   i_ms[1] = 38;
   j_ms[1] = 39;
   i_ms[2] = 39;
   j_ms[2] = 40;
   i_ms[3] = 39;
   j_ms[3] = 38;

   k_ms[0] =  0;   
   k_ms[1] =  0;
   k_ms[2] =  0;
   k_ms[3] =  0;

   for(int i = 0; i < num_ms; ++i) {
       x_ms[i] = x1[i_ms[i]];
       y_ms[i] = x2[j_ms[i]];
       z_ms[i] = x3[k_ms[i]];
   } 
}

void MasterPoint::display(){
   cout<<"MasterPoint::display(): " << endl;
   for(int k = 0; k < num_k_ms; k++) {
       for(int j = 0; j < num_j_ms; j++) {
           for(int i = 0; i < num_i_ms; i++) {
               int ijk = getIndex(i, j, k);
               cout << i_ms[ijk] << ' ' << j_ms[ijk] << ' '<<k_ms[ijk]<<' ';
               cout << x_ms[ijk] << ' ' << y_ms[ijk] << ' '<<z_ms[ijk] 
		    << endl;
	   }
       } 
   }          
}

void MasterPoint::output(){
   ofstream os("master_pts.out", ios::out);
// ofstream os("master_pts.out", ios::app);
   os << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(4);
   for(int k = 0; k < num_k_ms; k++) {
       for(int j = 0; j < num_j_ms; j++) {
           for(int i = 0; i < num_i_ms; i++) {
               int ijk = getIndex(i, j, k);
               os << i_ms[ijk] << ' ' << j_ms[ijk] << ' '<<k_ms[ijk]<<' ';
               os << x_ms[ijk] << ' ' << y_ms[ijk] << ' '<<z_ms[ijk] 
		  << endl;
	   }
       } 
   }          
   os<<endl;
   os.close();
}
