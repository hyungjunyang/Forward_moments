/* 
 *   double *Y, *dY[3]  -- Y=mean of ln(k/mu), dY[0]= dY/dx, dY[1]=dY/dy ...
 *                              size are nNode and 3*nNode
 *           Note: Y is changed to exp(-Y)=mu/k to avoid repeat calc..
 *   double *CY         -- covariance of Y, size is nNode^2 
 */

#ifndef _PERM_H
#define _PERM_H

#include <iostream>
#include <fstream>

#include "Common.h"
#include "Grid.h"
#include "Region.h"

#define XDIR 1
#define YDIR 2
#define ZDIR 3

class Perm {
	
public:
  Perm(Unit u, Region*, Grid*, char *Dir, bool data_cond, bool);
  virtual ~Perm();

  double getYAvg(int i) {return Y_Avg[i];}
  double getYVar(int i) {return Y_Var[i];}
  double getYAvg(int i, int j, int k) {
	         return Y_Avg[gridPtr->getPermIndex(i,j,k)];
                }
  double getYVar(int i, int j, int k) {
	         return Y_Var[gridPtr->getPermIndex(i,j,k)];
                }

  double getCorYY(int &i1, int &j1, int &k1, int &i2, int &j2, int &k2);
  double getCYY(int &i1, int &j1, int &k1, int &i2, int &j2, int &k2);
  double getCYY(int &ijk1, int &ijk2); // Pipatl
  void Perm_x( int i, int j, int k, double *permx );
  void Perm_y( int i, int j, int k, double *permy );
  void Perm_z( int i, int j, int k, double *permz );
  
  void PermCondition(int num_ms, double* weigh,
                     double* dY, double* dY_Std, int*);

  // set functions
  void readCondYMomntFile();
  void setCYY_Cond(int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double);
  void setCYY_Cond(int &ijk1, int &ijk2, double value);
  void setYAvg(int i, double value) { Y_Avg[i] = value;
	                              K_Avg[i] = exp(Y_Avg[i]);
                                    }
  void setYVar(int i, double value) { Y_Var[i] = value; }
  void setKAvg(int i) { K_Avg[i] = exp(Y_Avg[i]+ Y_Var[i]/2.); }
  void setZeroAtPts(int num_ms, double *weigh, double *dY_std);
  void setYVar(int i, int j, int k, double std) {
    Y_Var[gridPtr->getPermIndex(i,j,k)] = std * std;
  }
  
  void Copy_PermYAvgToTmp();
  void SetPermYAvg(int mst, int num_ms, double dY_avg, double *weigh);
  void Copy_TmpToPermYAvg();

  void wrtPermMomnt();
  void wrtPermMomnt(char *file);
private:
  bool debug;

  Unit unit;  
  char   *directory;
  Grid   *gridPtr;
  Region *regnPtr;

  // fields
  int i_ref, j_ref, k_ref;
  double *Y_Avg, *Y_Var, *K_Avg;
  double *Y_Var_Cond;
  double *Y_Avg_Tmp;

  void   initialize();
  void   addPermMomnt();
  void   delPermMomnt();
  void   setPermMomnt();
  double getExpoTrans(double YPerm);

  void   convertUnit();
  // temporal use;
  bool data_cond_;
  double * cyy_cond;
};  

#endif
