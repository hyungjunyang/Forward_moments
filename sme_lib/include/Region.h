/*
 * rInd : the index of regions
 *
 */
#ifndef _REGION_H
#define _REGION_H

#include "Control.h"
#include "Grid.h"
#include "Covariance.h"
#include "GaussCovariance.h"
#include "ExpCovariance.h"

class Region {
public:
  Region(ifstream &, Grid *gridPtr, Unit unit );
  ~Region(); 
  int   getNumOfRegion() {return nRegions; }

  // i - the index in the perm grid
  double getY_Avg(int& i) { return yMean[ rFineInd[i] ]; }
  double getY_Var(int& i) { return yVari[ rFineInd[i] ]; }
  double getCor(int &i1, int &j1, int &k1, int &i2, int &j2, int &k2);
 
  // i - the index in the pressure grid
  double getNodePermAvg(int i ) { return yMean[ rInd[i] ]; }
  double getNodePermVar(int i ) { return yVari[ rInd[i] ]; }
  double getNodePermCor(int i , int j ) { return yCorr[j][ rInd[i] ]; }
	 
private:
  Unit unit;
  Grid * gridPtr;

  double *x_Perm, *y_Perm, *z_Perm;
 
  int nRegions;
  int *rx[2], *ry[2], *rz[2];
  int    *yTrendType, *CovarianceType;
  double *yMean, *yVari, *yCorr[3], *azimuth;  
  double *yTrendInfo[9];
  ifstream *is_perm;
  int *nperm;
  Covariance **corPtr;
  int *rInd, *rFineInd;
  int typeRegions;
  ifstream *is_geology; 
  //interface
  int nInterfaces, *InterfaceDir, *InterfaceLoc;
  int *FirstRange[2], *SecondRange[2];
 
  void readDeck( ifstream & );
  void initialize();
  void assignIndex();
  void addCorObj();
  void initialize_interfaces();
  void convertUnit();
  void zeroingArray(int *array, int size);//(pipatl)To initialize new array(primarily for rFineInd)

  // debug part
  bool debug;
  void wrtRegionGridPermIndex();
};

#endif
