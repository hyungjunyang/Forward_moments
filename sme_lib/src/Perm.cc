#include "Perm.h"

Perm::Perm(Unit u, Region* reg, Grid* grid, char * Dir, bool data_cond, bool debug_ ) 
 : unit( u ), regnPtr ( reg ) , gridPtr ( grid ), 
   directory (Dir), debug(debug_) {
           
   data_cond_ =  data_cond;

   if(debug)  
      cout << endl << "Perm::Perm() " << endl;
  
   initialize();

   // cout << "Computing Y_Avg and Y_Var... " << endl; 
   setPermMomnt(); 
   
   if(unit == FIELD) convertUnit();
  
   if(debug) { 
      wrtPermMomnt();
   }
}

Perm::~Perm() {
   if(debug) 
      cout << "Perm::~Perm()" << endl;
   delPermMomnt();
}

void Perm::initialize() {
  i_ref = gridPtr->getPermNx()/2 - 1;
  j_ref = gridPtr->getPermNy()/2 - 1;
  k_ref = 0;
  addPermMomnt();
}

void Perm::addPermMomnt() {
  Y_Avg     = new double[ gridPtr->getNumPermNode() ];
  Y_Var     = new double[ gridPtr->getNumPermNode() ];
  K_Avg     = new double[ gridPtr->getNumPermNode() ];
  Y_Var_Cond= new double[ gridPtr->getNumPermNode() ];
  Y_Avg_Tmp = new double[ gridPtr->getNumPermNode() ];
  if( data_cond_ ) {
      int num_PermNode = gridPtr->getNumPermNode();
      cyy_cond  = new double[ num_PermNode * num_PermNode ];  
  } else {
      cyy_cond  = NULL;  
  }
}

void Perm::delPermMomnt() {
  delete[] Y_Avg;
  delete[] Y_Var;
  delete[] K_Avg;
  delete[] Y_Var_Cond;
  delete[] Y_Avg_Tmp;
  delete[] cyy_cond;
}

//===== setpermMomnt() =====
void Perm::setPermMomnt() {
  for(int k = 0; k < gridPtr->getPermNz(); ++k) {
     for(int j = 0; j < gridPtr->getPermNy(); ++j) {
         for(int i = 0; i < gridPtr->getPermNx(); ++i) {
             int ijk = gridPtr->getPermIndex(i, j, k);
             Y_Avg[ijk] = regnPtr->getY_Avg(ijk);
             Y_Var[ijk] = regnPtr->getY_Var(ijk);
             K_Avg[ijk] = getExpoTrans(Y_Avg[ijk]);
             //cout<<"NOTE:TRYING -- SPECIAL CAUTION!" << endl;
             //K_Avg[ijk] = getExpoTrans(Y_Avg[ijk] + 0.5 * Y_Var[ijk]);
             Y_Var_Cond[ijk] = Y_Var[ijk];
         } 
     }
  }
  if( data_cond_ ) {
      int num_PermNode = gridPtr->getNumPermNode();
      for(int k2 = 0; k2 < gridPtr->getPermNz(); ++k2) {
          for(int j2 = 0; j2 < gridPtr->getPermNy(); ++j2) {
              for(int i2 = 0; i2 < gridPtr->getPermNx(); ++i2) {
                  int ijk2 = gridPtr->getPermIndex(i2, j2, k2);
                  for(int k1 = 0; k1 < gridPtr->getPermNz(); ++k1) {
                      for(int j1 = 0; j1 < gridPtr->getPermNy(); ++j1) {
                          for(int i1 = 0; i1 < gridPtr->getPermNx(); ++i1) {
                              int ijk1 = gridPtr->getPermIndex(i1, j1, k1);
                              cyy_cond[ijk1 + ijk2 * num_PermNode]
                                  = sqrt(Y_Var[ijk1] * Y_Var[ijk2])
                                  * regnPtr->getCor(i1, j1, k1, i2, j2, k2); 
                          }
                      }          
                  }
              }
          }
      }
  } 
}

//===== ( readCondYMomntFile ) =====
void Perm::readCondYMomntFile() {
  cout << "READING CONDITIONAL PERMEABILITY MOMENTS DUE TO PERMEABILITY MEASUREMENTS" << endl;
  cout << "If no conditional permeability file, comment out 'readCondYMomntFile' from Flow.cc" << endl;
  //Set the system to be conditional
  data_cond_ = true;

  //Open files
  char fname[1000];
  ifstream is_ymn, is_cyy;
  sprintf(fname, "%sYmnCond.in", directory);
  is_ymn.open(fname, ios::in);
  if(!is_ymn.is_open()){
    cerr << "Cannot open " << fname << endl;
    exit(8);
  }

  sprintf(fname, "%sCYYCond.in", directory);
  is_cyy.open(fname, ios::in);
  if(!is_cyy.is_open()){
    cerr << "Cannot open " << fname << endl;
    exit(8);
  }

  //Read conditional permeability mean
  for(int k = 0; k < gridPtr->getPermNz(); ++k) {
    for(int j = 0; j < gridPtr->getPermNy(); ++j) {
      for(int i = 0; i < gridPtr->getPermNx(); ++i) {
        int ijk = gridPtr->getPermIndex(i, j, k);
        is_ymn >> Y_Avg[ijk];
        K_Avg[ijk] = getExpoTrans(Y_Avg[ijk]);
      }
    }
  }

  //Read C_YY (YVar)i
  int num_PermNode = gridPtr->getNumPermNode();
  for(int ijk1 = 0; ijk1 < num_PermNode; ++ijk1){
    for(int ijk2 = ijk1; ijk2 < num_PermNode; ++ijk2){
      is_cyy >> cyy_cond[ijk1 + ijk2 * num_PermNode];
      if(ijk1 != ijk2){
        cyy_cond[ijk2 + ijk1 * num_PermNode] = cyy_cond[ijk1 + ijk2 * num_PermNode];
      }
      else{
        Y_Var[ijk1] = cyy_cond[ijk1 + ijk2 * num_PermNode];
        Y_Var_Cond[ijk1] = Y_Var[ijk1];
      }
    }
  }

  is_ymn.close();
  is_cyy.close();
}

//===== ( getExpoTrans ) =====
double Perm::getExpoTrans(double YPerm) {
  return exp( YPerm );
}

//===== () =====
void Perm::Copy_PermYAvgToTmp(){
  for(int i = 0; i < gridPtr->getNumPermNode(); ++i) {
      Y_Avg_Tmp[i] = Y_Avg[i];
  }
}

void Perm::SetPermYAvg(int mst, int num_ms, double dY_avg, double *weigh ){
  for(int i = 0; i < gridPtr->getNumPermNode(); ++i) {
      Y_Avg[i] += weigh[ mst + num_ms * i] * dY_avg;
      K_Avg[i]  = getExpoTrans( Y_Avg[i] );
  }
}

void Perm::Copy_TmpToPermYAvg(){
  for(int i = 0; i < gridPtr->getNumPermNode(); ++i) {
      Y_Avg[i] = Y_Avg_Tmp[i];
      K_Avg[i] = getExpoTrans( Y_Avg[i] );
  }
}

// ===== PermCondition() =====
void Perm::PermCondition(int num_ms, double *weigh,
  double *dY_avg, double *dY_std, int* ijk_ms) {
  double Y_std_ijk;
  double sum_dY_Avg, sum_dY_Var;
  for(int k = 0; k < gridPtr->getPermNz(); ++k) {
     for(int j = 0; j < gridPtr->getPermNy(); ++j) {
         for(int i = 0; i < gridPtr->getPermNx(); ++i) {
             int ijk = gridPtr->getPermIndex(i, j, k);
             sum_dY_Avg = 0.;
             sum_dY_Var = 0.;
             for(int m = 0; m < num_ms; ++m) {
                sum_dY_Avg += weigh[ m + num_ms * ijk] * dY_avg[m];
                sum_dY_Var += weigh[ m + num_ms * ijk] * dY_std[m];
             }
             Y_Avg[ijk] += sum_dY_Avg;
             K_Avg[ijk]  = getExpoTrans(Y_Avg[ijk]);
             Y_std_ijk   = sqrt( Y_Var[ijk] );
             Y_std_ijk  +=  sum_dY_Var;             
             if(Y_std_ijk < 0) {
                Y_std_ijk = 0.;
             } else if( Y_std_ijk > 1.) {
                Y_std_ijk = 1.;
             }
             Y_Var[ijk] = Y_std_ijk * Y_std_ijk;
         } 
     }
  }
}

void Perm::setZeroAtPts(int num_ms, double *weigh, double *dY_std) {
  double Y_std_ijk, sum_dY_Var;
  for(int k = 0; k < gridPtr->getPermNz(); ++k) {
     for(int j = 0; j < gridPtr->getPermNy(); ++j) {
         for(int i = 0; i < gridPtr->getPermNx(); ++i) {
             int ijk = gridPtr->getPermIndex(i, j, k);
             sum_dY_Var = 0.;
             for(int m = 0; m < num_ms; ++m) {
                 sum_dY_Var += weigh[ m + num_ms * ijk] * dY_std[m];
             }
             Y_std_ijk   = sqrt( Y_Var[ijk] );
             Y_std_ijk  +=  sum_dY_Var;
             //if(Y_std_ijk < 0) {
             //   Y_std_ijk = 0.;
             //}
             //Y_Var[ijk] = Y_std_ijk * Y_std_ijk;
             Y_Var[ijk] = Y_std_ijk;
         } 
     }
  }
}

//===== (getCYY) =====
double Perm::getCYY(int &i1, int &j1, int &k1, int &i2, int &j2, int &k2) {
  int ijk1 = gridPtr->getPermIndex(i1, j1, k1);
  int ijk2 = gridPtr->getPermIndex(i2, j2, k2);
  double CYY_UnCd = sqrt(Y_Var[ijk1] * Y_Var[ijk2]) * regnPtr->getCor(i1, j1, k1, i2, j2, k2);
  if( data_cond_ ) {
      return cyy_cond[ijk1 + ijk2 * gridPtr->getNumPermNode()];
  } else {
      return CYY_UnCd;
  }
}

double Perm::getCYY(int &ijk1, int &ijk2)
{
	if(data_cond_){return cyy_cond[ijk1 + ijk2 * gridPtr->getNumPermNode()];}
	else 
	{
		cerr<<"Perm::getCYY(int &ijk1, int &ijk2) is for conditioning CYY only !!!";
	}

}

void Perm::setCYY_Cond(int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double value) {
  if( !data_cond_ ) {
      cout << "WRONG! in Perm::setCYY_Cond " << endl;
      exit(0);
  }
  int ijk1 = gridPtr->getPermIndex(i1, j1, k1);
  int ijk2 = gridPtr->getPermIndex(i2, j2, k2);
  cyy_cond[ijk1 + ijk2 * gridPtr->getNumPermNode()] -= value;
}

void Perm::setCYY_Cond(int &ijk1, int &ijk2, double value) {
  if( !data_cond_ ) {
      cout << "WRONG! in Perm::setCYY_Cond " << endl;
      exit(0);
  }
  cyy_cond[ijk1 + ijk2 * gridPtr->getNumPermNode()] -= value;
}


//===== (getCorYY) =====
double Perm::getCorYY(int &i1, int &j1, int &k1, int &i2, int &j2, int &k2) {
  return regnPtr->getCor(i1, j1, k1, i2, j2, k2);
}

//===== get Perm_x (permeability at the x - direction) =====
void Perm::Perm_x( int i, int j, int k, double *perm_x ) {
  int i_Perm      = gridPtr->toIPerm(i);
  int j_Perm      = gridPtr->toJPerm(j);
  int k_Perm      = gridPtr->toKPerm(k);
  int ijk_Perm_n1 = gridPtr->getPermIndex(i_Perm - 1, j_Perm, k_Perm);
  int ijk_Perm    = gridPtr->getPermIndex(i_Perm    , j_Perm, k_Perm);
  int ijk_Perm_p1 = gridPtr->getPermIndex(i_Perm + 1, j_Perm, k_Perm);
  perm_x[0]       = K_Avg[ijk_Perm_n1];
  perm_x[1]       = K_Avg[ijk_Perm   ];
  perm_x[2]       = K_Avg[ijk_Perm_p1];
}

//===== get Perm_y (permeability at the y - direction) =====
void Perm::Perm_y( int i, int j, int k, double *perm_y ) {
  int i_Perm      = gridPtr->toIPerm(i);
  int j_Perm      = gridPtr->toJPerm(j);
  int k_Perm      = gridPtr->toKPerm(k);
  int ijk_Perm_n1 = gridPtr->getPermIndex(i_Perm, j_Perm - 1, k_Perm);
  int ijk_Perm    = gridPtr->getPermIndex(i_Perm, j_Perm    , k_Perm);
  int ijk_Perm_p1 = gridPtr->getPermIndex(i_Perm, j_Perm + 1, k_Perm);
  perm_y[0]       = K_Avg[ijk_Perm_n1];
  perm_y[1]       = K_Avg[ijk_Perm   ];
  perm_y[2]       = K_Avg[ijk_Perm_p1];
}

//===== get Perm_z (permeability at the z - direction) =====
void Perm::Perm_z( int i, int j, int k, double *perm_z ) {
  int i_Perm      = gridPtr->toIPerm(i);
  int j_Perm      = gridPtr->toJPerm(j);
  int k_Perm      = gridPtr->toKPerm(k);
  int ijk_Perm_n1 = gridPtr->getPermIndex(i_Perm, j_Perm, k_Perm - 1);
  int ijk_Perm    = gridPtr->getPermIndex(i_Perm, j_Perm, k_Perm    );
  int ijk_Perm_p1 = gridPtr->getPermIndex(i_Perm, j_Perm, k_Perm + 1);
  perm_z[0]       = K_Avg[ijk_Perm_n1];
  perm_z[1]       = K_Avg[ijk_Perm   ];
  perm_z[2]       = K_Avg[ijk_Perm_p1];
}


//===== ( writeCYY ) =====
void Perm::wrtPermMomnt(){
  char outfile1[80], outfile2[80];
  sprintf( outfile1, "%sCYY_1D_XY.out", directory );
  sprintf( outfile2, "%sCYY_2D_XY.out", directory );
  //cout << outfile1 << endl;
  //cout << outfile2 << endl;
  ofstream os1(outfile1, ios::out);
  ofstream os2(outfile2, ios::out);
  
  os1 << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(4);
  os2 << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(4);
 
  int ijk;
  if(gridPtr->getNx() > 1 && gridPtr->getNy() > 1) {
     int k = k_ref;
     for(int j = 0; j < gridPtr->getPermNy(); j++) {
         int i = j;    
         ijk = gridPtr->getPermIndex(i, j, k);
         os1<< i << ' ' 
            << j << ' '
            << Y_Avg[ijk]<<' '
            << K_Avg[ijk]<<' '
            << Y_Var[ijk]<<' '
            << getCYY(i, j, k, i_ref, j_ref, k_ref) << endl;
     }
  }
  
  if(gridPtr->getNx() > 1 && gridPtr->getNy() > 1) {
     int k = k_ref;
     for(int j = 0; j < gridPtr->getPermNy(); j++) {
         for(int i = 0; i < gridPtr->getPermNx(); i++) {
             ijk = gridPtr->getPermIndex(i, j, k);
             os2<< i << ' ' 
                << j << ' '
                << Y_Avg[ijk]<<' '
                << K_Avg[ijk]<<' '
                << Y_Var[ijk]<<' '
                << getCYY(i, j, k, i_ref, j_ref, k_ref) << endl;
         } 
         os2 << endl;
     }
  }
  os1.close();
  os2.close();
}

void Perm::wrtPermMomnt(char *file) {
  char outfile1[80], outfile2[80], outfile3[80];
  sprintf( outfile1, "%sCYY_1D_XY.out", file );
  sprintf( outfile2, "%sCYY_2D_XY.out", file );
  sprintf( outfile3, "%sCYY_Mat_XY.out", file);
  ofstream os1(outfile1, ios::out);
  ofstream os2(outfile2, ios::out);
  ofstream os3(outfile3, ios::out);
  
  os1 << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(4);
  os2 << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(4);
  os3 << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(4);
 
  int ijk;
  if(gridPtr->getNx() > 1 && gridPtr->getNy() > 1) {
     int k = k_ref;
     for(int j = 0; j < gridPtr->getPermNy(); j++) {
         int i = j;    
         ijk = gridPtr->getPermIndex(i, j, k);
         os1<< i << ' ' 
            << j << ' '
            << Y_Avg[ijk]<<' '
            << K_Avg[ijk]<<' '
            << getCYY(i, j, k, i, j, k) << ' '
            << getCYY(i, j, k, i_ref, j_ref, k_ref) << endl;
     }
  }
  
  if(gridPtr->getNx() > 1 && gridPtr->getNy() > 1) {
     int k = k_ref;
     for(int j = 0; j < gridPtr->getPermNy(); j++) {
         for(int i = 0; i < gridPtr->getPermNx(); i++) {
             ijk = gridPtr->getPermIndex(i, j, k);
             os2<< i << ' ' 
                << j << ' '
                << Y_Avg[ijk]<<' '
                << K_Avg[ijk]<<' '
                << getCYY(i, j, k, i, j, k) << ' '
                << getCYY(i, j, k, i_ref, j_ref, k_ref) << endl;
         } 
         os2 << endl;
     }
  }

  if(gridPtr->getNx() > 1 && gridPtr->getNy() > 1) {
     //2D full CYY matrix
     for(int j2 = 0; j2 < gridPtr->getPermNy(); j2++) {
	for(int i2 = 0; i2 < gridPtr->getPermNx(); i2++) {
           for(int j1 = 0; j1 < gridPtr->getPermNy(); j1++) {
              for(int i1 = 0; i1 < gridPtr->getPermNx(); i1++) {
                 os3<< getCYY(i1, j1, k_ref, i2, j2, k_ref) << ' ';
              }
           }
           os3 << endl;
        }
     }
  }
  os1.close();
  os2.close();
  os3.close();
}

//===== ( convertUnit ) =====
// NOTE: IMPORTANT
// This conversion is different from the original one
// the original Y = mu/k
// Here K/mu        
//how do you convert C_YY
void Perm::convertUnit() {
  double aaa = md_m2 / cp_pasec;
  for(int i = 0; i < gridPtr->getNumPermNode(); ++i) {
      Y_Avg[i] += log(md_m2);
      K_Avg[i] *= aaa;
      //cout <<Y_Avg[i] << ' ' << K_Avg[i]  << endl;
  }
  //exit(0);
}
