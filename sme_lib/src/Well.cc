/*
 * File: Well.cc
 * 
 * --------------------------------------------------------------------
 *  WConType wConType -- well control type (1:RATE; 2:PRES) 
 *  int wLens         -- number of penetrations
 *  int *wLoc[3]      -- i,j,k index for each well block
 *  double *wTrans    -- well transmissiblity for each well block
 *  double wStr       -- pres for const pres well, rate for const rate well
 *  
 * --------------------------------------------------------------------
 * Implementation for Well class
 */

#include "Well.h"

/* --- Public Methods ------------------------------ */
// --- Constructors and destructor ------------
Well::Well() {
  debug = true;
  //debug = false;
  if(debug) cout << "Well::Well()" << endl;
}

Well::~Well() {
  if(debug) cout << "Well::~Well()" << endl;
  for(int i = 0; i < 3; i++) delete[] wLocIJK[i];
  for(int i = 0; i < 3; i++) delete[] wLocXYZ[i];
  delete[] wTrans; 
}

void Well::readData(ifstream &is, Unit unit) {
  int i;
  // --- well control type and well strength --
  Junk(is); is >> i >> wStr;
  if(i==1){wStr *= -1;} //To make production rate positive
  Junk(is); is >> isInjcWell >> isProdWell >> angle[0] >> angle[1] >> radius >> npts;
  if(debug) {
     cout << isInjcWell <<' '<< isProdWell <<' '<< angle[0] << ' ';
     cout << angle[1]   <<' '<< radius     <<' '<< npts << endl;
  }
  wConType = (WConType) i;

  // --- number of well blocks (number of penetration ) ---
  Junk(is); is >> wLens;
  if(debug) cout << "wLens = " << wLens << endl;
  wTrans = new double[wLens];
  for(i = 0; i < 3; i++) {
      wLocIJK[i] = new int   [wLens];
      wLocXYZ[i] = new double[wLens];
  } 
  // --- index and transmissibility of well blocks ---
  for(i = 0; i < wLens; i++) {
      Junk(is); 
      is >> wLocIJK[0][i] >> wLocIJK[1][i] >> wLocIJK[2][i] >> wTrans[i];
      if(debug) cout <<            i  << ' '
	             << wLocIJK[0][i] << ' '
                     << wLocIJK[1][i] << ' '
                     << wLocIJK[2][i] << ' '
                     << wTrans[i] << ' '
                     << endl;
  }
  for(i = 0; i < wLens; i++) {
      Junk(is);
      is >> wLocXYZ[0][i] >> wLocXYZ[1][i] >> wLocXYZ[2][i];
      if(debug) cout <<            i  << ' '
	             << wLocXYZ[0][i] << ' '
                     << wLocXYZ[1][i] << ' '
                     << wLocXYZ[2][i] << endl;
  }
  // --- return back to "0" based index ------
  for(i = 0; i < wLens; i++) { 
      wLocIJK[0][i]--; 
      wLocIJK[1][i]--; 
      wLocIJK[2][i]--; 
  }

  if(unit==FIELD) {
     convertUnit();
     if(wConType == RATE ) {
	cout << "wStr of rate = " << wStr << endl;
     }
  }
}

// --- output the well  -------------------------------------
ostream &operator<< (ostream &os, Well &w)
{ 
  (w.wConType==RATE)? os << "Fixed Flow Rate = " : os << "Fixed Pressure = ";
  os << w.wStr << endl; 

  for(int i=0; i<w.wLens; i++) { 
      os << "At Loc Index = ("
         << w.wLocIJK[0][i]<<", "
         << w.wLocIJK[1][i]<<", "
         << w.wLocIJK[2][i]<<"), or " << endl; 
      os << "Physical Loc = ("
         << w.wLocXYZ[0][i]<<", "
         << w.wLocXYZ[1][i]<<", "
         << w.wLocXYZ[2][i]<<")"<<endl; 
      os << "Trans = " << w.wTrans[i]<<endl;
  }
  return os;
}

/* --- Private Methods ----------------------------- */
// --- convert field unit to metric unit -----------
void Well::convertUnit() {
  wStr *= (wConType == RATE)? stbPday_m3Psec : psia_pa;

  double aaa = stbPday_m3Psec/psia_pa;
  //cout << "aaa =  "<< aaa << endl;
  for(int i = 0; i < wLens; i++) {
      //cout << "wTrans[i]  = " << wTrans[i] << endl;	  
      wTrans[i] *= aaa;
      //cout << "wTrans[i]  = " << wTrans[i] << endl;
  }
  for(int i = 0; i < wLens; i++) {
      wLocXYZ[0][i] *= ft_m;
      wLocXYZ[1][i] *= ft_m;
      wLocXYZ[2][i] *= ft_m;
  }
  radius *= ft_m;
}
