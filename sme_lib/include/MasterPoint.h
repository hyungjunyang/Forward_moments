/*c   seed     the seed number, readed from the input file             !
c     num_i_ms no. of columms to be generated for master points        !
c     num_j_ms no. of columms to be generated for master points        !
c     nx,ny    the number of nodes in the x-, y- directions of perm    !
c     x1,x2    the coordinates of x and y directions of perm filed     !
c     j_ms     the j-th grid of the selected master point              !
c     x_ms     the x-coordinate of the master points                   !
c     y_ms     the y-coordinate of the master points                   !
c     z_ms     the z-coordinate of the master points                   !
*/
#ifndef _MASTERPOINT_H
#define _MASTERPOINT_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <math.h>

float ran3(long *idum);

class MasterPoint{
public:
  MasterPoint(long, int num_i, int num_j, int num_k);
 ~MasterPoint();

  void initialization();
  int  getLength() {return num_ms;}
  int  getIndex(int& i, int& j, int& k) { return i + j * num_i_ms 
	                                           + k * num_i_ms * num_j_ms; }
  int* getIm()      { return i_ms; }
  int* getJm()      { return j_ms; }
  int* getKm()      { return k_ms; }
  int* getIJK()     { return ijk_ms;}

  int  getIm(int i) { return i_ms[i];}
  int  getJm(int i) { return j_ms[i];}
  int  getKm(int i) { return k_ms[i];}
  int  getIJK(int i) { return ijk_ms[i];}

  void setMstPts(double*, double*, double*);
  void calMasterPts(int nx, int ny, int nz, double* x1, double* x2, double* x3);
  void display();
  void output(); 

private:
  bool debug;
  long seed;
  int num_i_ms, num_j_ms, num_k_ms;
  int num_ms;
  int    *i_ms, *j_ms, *k_ms, *ijk_ms;
  double *x_ms, *y_ms, *z_ms; 

};
#endif
