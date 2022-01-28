/*
 * File: Control.cc
 * ----------------------------------
 * Implementation for Control class
 */
#include "Control.h"

Control::Control(ifstream &is, Grid* g, Unit unit)
: gridPtr(g) {

  debug = false;
  if(debug) cout << "Control::Control()" << endl;
  
  bType = new BType[6];
  bCond = new double[6];
  bLength = new double[6]; 
  bProd = new double[6]; 
  bInjc = new double[6]; 

  readData(is);
  if(unit == FIELD) convertUnit();
  calcBoundPoint(); 
 
}

Control::~Control() {
  if(debug) cout << "Control::~Control()" << endl;	
  delete[] dt;
  delete[] bType;
  delete[] bCond;
  delete[] bLength; 
  delete[] bPtSetXn; 
  delete[] bPtSetXp; 
  delete[] bPtSetYn;
  delete[] bPtSetYp; 
  delete[] bProd;
  delete[] bInjc;
}

//==========================================
ostream &operator<< (ostream &os, Control &c)
{
  int i;
  os << endl;
  os << "===== Boundary Conditions =====" << endl;
  os << "X+   X-   Y+   Y-   Z+   Z-" << endl;
  for (i = 0; i < 6; i++) os << c.bType[i] << "    ";  os << endl;
  for (i = 0; i < 6; i++) os << c.bCond[i] << "    ";  os << endl;
  
  os << endl;
  os << "===== Initial Pressure =====" << endl;
  os << "pInit = " << c.pInit << endl;
  
  os << endl;
  os << "===== Solver Info =====" << endl;
  os << "Solver Choice = " << c.solverChoice << endl;
 
  os << endl; 
  os << "===== Reference Point For Output =====" << endl;
  os << "(" << c.icc << ", " 
	    << c.jcc << ", " 
	    << c.kcc << ")\n";

  os << endl;
  os << "===== Timesteps Info =====" << endl;

  // transient problem
  if(c.dt[0] > 0.0) {
      os << "Total Number of Timesteps: " << c.nDt << endl;
      for (i = 0; i <c.nDt; i++) 
	   os << c.dt[i] << " ";  
      os << endl;
  }
  else {
      os << "Steady State Run, No Timesteps.\n";
  }
  os << endl;
  os << "Flag for Conditioning = " << c.flag_press_condition ;
  if(c.flag_press_condition == 0) {
     os<<" => No Conditioning" << endl;
  }else {
     os<<" => Using Pressure Measurements for Conditioning" << endl;
  }

  os << "Swc = " << c.Swc << ' '
     << "Sor = " << c.Sor << ' '
     << "viscos_w = " << c.viscos_w <<' ' 
     << "viscos_o = " << c.viscos_o <<' '
     << "Viscosity Radio = " << c.viscosR << endl ;

  os << "Time for Production = " << c.time_production <<' '
     << "Number of time steps = " << c.num_time_step << endl;
  
  return os;
}
//==== caculate the starting points of SL at boundary 
//==== now only support 2D case 
void Control::calcBoundPoint() {

    int nx = gridPtr->getNx();
    int ny = gridPtr->getNy(); 
    double* bdx = gridPtr->getBdx(); 
    double* bdy = gridPtr->getBdy(); 
    double* x = gridPtr->getX(); 
    double* y = gridPtr->getY(); 
    double Lx = gridPtr->getDim( 0 ); 
    double Ly = gridPtr->getDim( 1 ); 

    int i_block = 0; 
    int j_block = 0;
    double x_pt = 0.0;
    double y_pt = 0.0; 
    double del_x = 0.0; 
    double del_y = 0.0; 

    // x-direction boundary 
    if ( bLength[0] > 0 ) { 
       int len = bLength[0]; 
       bPtSetXn = new Point[ len ]; 
       i_block = 0; 
       x_pt = 0.0; 
       del_y = Ly / bLength[0]; 
       for ( int j = 0; j < bLength[0]; j++ ) { 
          y_pt = ( j + 0.5 ) * del_y;  
          for ( int j_tmp = 0; j_tmp < ny; j_tmp++ ) {
              if ( y[j_tmp] - 0.5 * bdy[j_tmp] <= y_pt && y_pt < y[j_tmp] + 0.5 * bdy[j_tmp] ){
                  j_block = j_tmp;    
              } 
          }
          Point tmp_pt( x_pt, y_pt, i_block, j_block ); 
          bPtSetXn[j] = tmp_pt;
       }
    }
    if ( bLength[1] > 0 ) {
       int len = bLength[1];  
       bPtSetXp = new Point[ len ]; 
       i_block = nx - 1; 
       x_pt = Lx;
       del_y = Ly / bLength[1];
       for ( int j = 0; j < bLength[1]; j++ ) {
          y_pt = ( j + 0.5 ) * del_y;
          for ( int j_tmp = 0; j_tmp < ny; j_tmp++ ) {
              if ( y[j_tmp] - 0.5 * bdy[j_tmp] <= y_pt && y_pt < y[j_tmp] + 0.5 * bdy[j_tmp] ){
                  j_block = j_tmp;
              }
          }
          Point tmp_pt( x_pt, y_pt, i_block, j_block );
          bPtSetXp[j] = tmp_pt;
       }
    }
    // y-direction boundary 
    if ( bLength[2] > 0 ) { 
       int len = bLength[2]; 
       bPtSetYn = new Point[ len ]; 
       j_block = 0; 
       y_pt = 0.0; 
       del_x = Lx / bLength[2]; 
       for ( int i = 0; i < bLength[2]; i++ ) { 
          x_pt = ( i + 0.5 ) * del_x;  
          for ( int i_tmp = 0; i_tmp < nx; i_tmp++ ) {
              if ( x[i_tmp] - 0.5 * bdx[i_tmp] <= x_pt && x_pt < x[i_tmp] + 0.5 * bdx[i_tmp] ){
                  i_block = i_tmp;    
              } 
          }
          Point tmp_pt( x_pt, y_pt, i_block, j_block ); 
          bPtSetYn[i] = tmp_pt;
       }
    }
    if ( bLength[3] > 0 ) { 
       int len = bLength[3]; 
       bPtSetYp = new Point[ len ]; 
       j_block = ny - 1; 
       y_pt = Ly; 
       del_x = Lx / bLength[3]; 
       for ( int i = 0; i < bLength[3]; i++ ) { 
          x_pt = ( i + 0.5 ) * del_x;  
          for ( int i_tmp = 0; i_tmp < nx; i_tmp++ ) {
              if ( x[i_tmp] - 0.5 * bdx[i_tmp] <= x_pt && x_pt < x[i_tmp] + 0.5 * bdx[i_tmp] ){
                  i_block = i_tmp;    
              } 
          }
          Point tmp_pt( x_pt, y_pt, i_block, j_block ); 
          bPtSetYp[i] = tmp_pt;
       }
    }

}
//===== (getBoundary points) ===
Point* Control::getBoundPts( int i ){ 
    switch (i) {
        case 0: 
            return bPtSetXn; 
            break; 
        case 1: 
            return bPtSetXp; 
            break; 
        case 2:
            return bPtSetYn; 
            break; 
        case 3:
            return bPtSetYp; 
    }
    return NULL; 
}

//===== (readData) =====
void Control::readData(ifstream &is) {
  int i, j;
  // --- boundary conditions --- 
  for(i = 0; i < 6; i++) {
      Junk(is); is >> j >> bCond[i] >> bInjc[i] >> bProd[i] >> bLength[i];
      bType[i] = (BType) j;
  }

  // --- initial condition ---
  Junk(is); is >> pInit;
  
  // --- solver choice ---
  Junk(is); is >> j;  solverChoice = (SolverType) j;
  if(debug) cerr << "Solver type: " << solverChoice << endl;   

  // --- reference point for output ---
  Junk(is); is >> icc >> jcc >> kcc;
  icc--; jcc--; kcc--;   // back to 0 based index

  // --- timestep control ---
  Junk(is); is >> nDt;
  dt = new double[nDt];
  int nDt_i;  double dt_i; i = 0;

  // read and assign the dt[]
  while(i < nDt) {
      Junk(is);  
      is >> nDt_i >> dt_i;
      for(j = 0; j < nDt_i; j++) 
	  dt[i+j] = dt_i;
      i += j;
  }
  
  // flag for pressure conditioning
  Junk(is); 
  is >> flag_press_condition;

  // flag for saturation conditioning (Pipatl)
  Junk(is);
  is >> flag_sat_condition;

  // flag for production conditioning (Pipatl)
  Junk(is);
  is >> flag_prod_condition;
 
  // production
  Junk(is);
  is >> Swc >> Sor >> viscos_w >> viscos_o;
  viscosR = viscos_o / viscos_w;

  Junk(is);
  is >> time_production >> num_time_step;
}

//===== convert the field unit to metric unit =========================
void Control::convertUnit() {
// cout << pInit << endl;
  pInit *= psia_pa;
// cout << pInit << endl;
// exit(0);

  for(int i = 0; i < nDt; i++) 
      dt[i] *= day_sec;
  for(int i=0; i<6; i++)
      bCond[i] *= (bType[i]==CONST_PRES)? psia_pa : stbPday_m3Psec;
}

