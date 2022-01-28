#include "Fluid.h"

//===== ( Fluid() ) =====
Fluid::Fluid(ifstream &is, Grid* g, Unit u) 
: unit(u), gridPtr(g) {

  debug = false;  
  if(debug) cout << "Fluid::Fluid()" << endl;
  
  readData( is );     // all initialization are
  initialize();

  if(unit == FIELD) convertUnit();
  update(NULL, true);
}

//===== ( ~Fluid() ) =====
Fluid::~Fluid() { 
  //cout<<"Fluid::~Fluid()"<<endl;
  delete[] poro;
  delete[] den_w;
  delete[] den_o;
  //delete[] a;
  //delete[] b;
  //delete[] c;
  delete gamma;
  delete den_cf;
  delete[] poro_crf;
}

//===== ( readData(ifstream &is) ) =====
void Fluid::readData(ifstream &is) {        
  cout << endl;
  cout << "===== Fluid Info =====" << endl; 
  Junk(is);
  is >> rPoro >> cr >> cf_w >> cf_o >> rDen_w >> rDen_o >> rPres;
  cout << "rPoro = " << rPoro << endl;
  cout << "cr, cf_w, cf_o = " 
       <<  cr   << ' '
       <<  cf_w << ' ' 
       <<  cf_o << endl;
  cout <<  "rDen_w, rDen_o = " 
       <<   rDen_w <<' '
       <<   rDen_o << endl;
  // --- Saturation profile (front) shape ----
  Junk(is);
  is >> satPType >> a0 >> a1; 
  cout << "satPType = " <<  satPType << endl;
  cout << "a0 = " << a0 << " a1 = " << a1<<endl;
}

//===== ( initialize() ) =====
void Fluid::initialize() {
  int nz = gridPtr->getNz();
  poro  =             new double[ gridPtr->getnNode() ];
  den_w = ( nz > 1) ? new double[ gridPtr->getnNode() ] : NULL;
  den_o = ( nz > 1) ? new double[ gridPtr->getnNode() ] : NULL;
  //a     = ( nz > 1) ? new double[ gridPtr->getnNode() ] : NULL;
  //b     = ( nz > 1) ? new double[ gridPtr->getnNode() ] : NULL;
  //c     =             new double[ gridPtr->getnNode() ];
  gamma    = ( nz > 1) ? new double[ gridPtr->getnNode() ] : NULL;
  den_cf   = ( nz > 1) ? new double[ gridPtr->getnNode() ] : NULL;
  poro_crf =             new double[ gridPtr->getnNode() ];
}

//===== ( Update ) =====
void Fluid::update(const double *P, bool init) {
  // --- calculate porosity and c ------------
  if( cr > 0 || init ) {
      for(int i = 0; i < gridPtr->getnNode(); i++) {
          poro[i]= (init)? rPoro : rPoro*exp( cr*(P[i]-rPres) );
          //c[i] = poro[i]*(cr+cf_o+cf_w);
          poro_crf[i] = poro[i]*(cr + cf_o + cf_w);          
      }
  }

  // --- calculate den_w, den_o, and a, b ----
  if( gridPtr->getNz() > 1 ) {
      if( cf_w > 0 || init ) {
          for(int i = 0; i < gridPtr->getnNode(); i++) {
              den_w[i] = (init) ? 
                         rDen_w : rDen_w * exp( cf_w*(P[i]-rPres) );
          }
      }
      if( cf_o>0 || init ) {
          for(int i = 0; i < gridPtr->getnNode(); i++) {
              den_o[i] = (init) ? 
                         rDen_o : rDen_o * exp( cf_o*(P[i]-rPres) );
          }
      }
      // add option for steady state problem, i.e., 
      // domain is fully saturated with water
      if( init ) {
          if( cf_w>0.0 || cf_o>0.0 || cr>0.0 ) 
            (satPType == 1) ? 
             usePhaseInfo_linear() : usePhaseInfo_eclipse();
          else {
             for(int k = 0; k < gridPtr->getNz(); k++) {
                 for(int j = 0; j < gridPtr->getNy(); j++) {
                     for(int i = 0; i < gridPtr->getNx(); i++) {
                         int ijk = gridPtr->getIndex( i, j, k );    
                         //a[ijk] = 9.81 * den_w[ijk];
                         //b[ijk] = 0.0;
                         //c[ijk] = 0.0;
                         gamma[ijk] = 9.81 * den_w[ijk];
                         den_cf[ijk] = 0.0;
                         poro_crf[ijk] = 0.0;
                     }
                 }
             }
          }
      }
  }
}

//===== ( usePhaseInfo_linear ) =====
// xx, zz are scaled coordinate, within [0, 1]
void Fluid::usePhaseInfo_linear() {
  double xx, zz, zl, krw, kro, sw;
  for(int k = 0; k < gridPtr->getNz(); k++) {  
      for(int j = 0; j < gridPtr->getNy(); j++) {
          for(int i = 0; i < gridPtr->getNx(); i++) {
              int ijk = gridPtr->getIndex( i, j, k );
              xx = gridPtr->getX(i) / gridPtr->getDim(0);
              zz = gridPtr->getZ(k) / gridPtr->getDim(2);
              zl = a0 + xx * a1;      // Sw profile's z correspoding to xx
              if(zz > zl) {           // ahead the front, the oil zone
                 //a[ijk] = 9.8 * den_o[ijk];
                 //b[ijk] = 2.0 * gamma[ijk] * cf_o;
                 //c[ijk] = poro[ijk] * ( cr + cf_o );
                 gamma[ijk] = 9.81 * den_w[ijk];
                 den_cf[ijk] = 2.0 * gamma[ijk] * cf_o;
                 poro_crf[ijk] = poro[ijk] * ( cr + cf_o );
              } else {
                 // behind the front, w+o
                 // distribute sw linearly in xx from xx=0 to the front
                 // --- at xx=0, sw=1; at front, sw=0.7 ---
                 sw = (a0 != zz)? 1-(a0-zl)/(a0-zz)*0.3 : 1.0;
                 krw = sw * sw;      // relative perms: krw and kro      
                 kro = ( 1 - sw ) * ( 1 - sw );
                 //a[ijk] = 9.8*(den_w[ijk]*krw+den_o[ijk]*kro);
                 //b[ijk] = 2.0*9.8*(den_w[ijk]*krw*cf_w+den_o[ijk]*kro*cf_o);
                 //c[ijk] = poro[ijk]*(cr+cf_o*(1-sw)+cf_w*sw);
                 gamma[ijk] = 9.8*(den_w[ijk]*krw+den_o[ijk]*kro);
                 den_cf[ijk] = 2.0*9.8*(den_w[ijk]*krw*cf_w+den_o[ijk]*kro*cf_o);
                 poro_crf[ijk] = poro[ijk]*(cr+cf_o*(1-sw)+cf_w*sw);
              }
          }
      }
  }
}

//===== ( usePhaseInfo_eclipse ) =====
void Fluid::usePhaseInfo_eclipse() {
        
  // xx, zz are scaled coordinate, within [0, 1]
  double xx, zz, rr, krw, kro, sw;
  for(int k = 0; k < gridPtr->getNz(); k++) {  
      for(int j = 0; j < gridPtr->getNy(); j++) {
          for(int i = 0; i < gridPtr->getNx(); i++) {
              int ijk = gridPtr->getIndex( i, j, k );
              xx  = gridPtr->getX(i) / gridPtr->getDim(0);                    
              zz  = gridPtr->getZ(k) / gridPtr->getDim(2);
              rr = sqrt(zz*zz/a0/a0+xx*xx/a1/a1);  
              // eclipse radus 
              // outside the eclipse(W front), oil zone
              if(rr > 1.0) {                         
                 //a[ijk] = 9.8 * den_o[ijk];
                 //b[ijk] = 2.0 * gamma[ijk] * cf_o;
                 //c[ijk] = poro[ijk] * ( cr + cf_o );
                 gamma[ijk] = 9.8 * den_o[ijk];
                 den_cf[ijk] = 2.0 * gamma[ijk] * cf_o;
                 poro_crf[ijk] = poro[ijk] * ( cr + cf_o );
              } else {                             
                 // inside the eclipse, w+o
                 // distribute sw linearly in rr from rr=0 to the front
                 // ar rr=0, sw=1; at front (rr=1) sw=0.7 ---
                 sw  = 1 - rr * 0.3;
                 krw = sw * sw;   // relative perms: krw and kro    
                 kro = (1 - sw) * (1 - sw);
                 //a[ijk] =     9.8*(  den_w[ijk] * krw 
                 //                  + den_o[ijk] * kro );
                 //b[ijk] = 2.0*9.8*(  den_w[ijk] * krw * cf_w
                 //                  + den_o[ijk] * kro * cf_o );
                 //c[ijk] = poro[ijk] * ( cr + cf_o * (1 - sw ) 
                 //                          + cf_w *      sw ); 
                 gamma[ijk] = 9.8*(  den_w[ijk] * krw 
                                   + den_o[ijk] * kro );
                 den_cf[ijk] = 2.0*9.8*(  den_w[ijk] * krw * cf_w
                                        + den_o[ijk] * kro * cf_o );
                 poro_crf[ijk] = poro[ijk] * (  cr + cf_o * (1 - sw ) 
                                              + cf_w * sw );
              }
          }
      }
  }
}

//===== ( convertUnit() ) =====
void Fluid::convertUnit() {
  cr     /= psia_pa; 
  cf_w   /= psia_pa; 
  cf_o   /= psia_pa; 
  rPres  *= psia_pa; 
  rDen_w *= lbmPft3_kgPm3;  
  rDen_o *= lbmPft3_kgPm3;
}

