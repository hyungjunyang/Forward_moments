/*
 * File: Grid.h
 * --------------------------------------------------------------------
 * This class defines the grid (point distributed or cell centered)
 *
 * --------------------------------------------------------------------
 **Geometry Data*******************************************************
 *   GridType gridType     -- (1:POINT_DIS; 2:CELL_CEN)
 *   int nx,ny,nz,nNode    -- x, y, z dimension, nNode=nx*ny*nz
 *   double *x, *y, *z     -- node coordinate, size are nx, ny, nz
 *   double *dx,*dy,*dz    -- distance between nodes, size are nx-1,ny-1,nz-1
 *                            dx[i]=x[i]-x[i-1], dy[j]=y[j]-y[j-1], ...
 *   double *bdx,*bdy,*bdz -- block(control volume) size, size are nx,ny,nz
 *   double dim[3]         -- physical dimension of the reservoir
 *
 * --------------------------------------------------------------------
 */

#ifndef _Grid_H
#define _Grid_H

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <math.h> 

#include "Common.h"

class Grid {
 public:
  // --- Constructor and destructor ---
   Grid( ifstream &input, Unit unit, char *Dir );
  ~Grid();

  int getNx() const { return nx; }
  int getNy() const { return ny; }
  int getNz() const { return nz; }
  int getnNode()    { return nNode; }
  
  double  getX(int i) const { return x[i]; }
  double  getY(int i) const { return y[i]; }
  double  getZ(int i) const { return z[i]; }
  double* getX(     ) const { return x; }
  double* getY(     ) const { return y; }
  double* getZ(     ) const { return z; }

  double* getDx(     ) const { return dx; }
  double* getDy(     ) const { return dy; }
  double* getDz(     ) const { return dz; }
  double  getDx(int i) const { return dx[i]; }
  double  getDy(int i) const { return dy[i]; }
  double  getDz(int i) const { return dz[i]; }
 
  double* getBdx(     ) const { return bdx; }
  double* getBdy(     ) const { return bdy; }
  double* getBdz(     ) const { return bdz; }
  double  getBdx(int i) const { return bdx[i]; }
  double  getBdy(int i) const { return bdy[i]; }
  double  getBdz(int i) const { return bdz[i]; }

  double getDim(int i) const {return dim[i];}

  
  int getPermNx() const { return nxPerm; }
  int getPermNy() const { return nyPerm; }
  int getPermNz() const { return nzPerm; }
  int getNumPermNode()  { return nNodePerm; }

  double* getPermX() const { return xPerm; }
  double* getPermY() const { return yPerm; }
  double* getPermZ() const { return zPerm; }
  double  getPermX(int& i)  { return xPerm[i]; }
  double  getPermY(int& i)  { return yPerm[i]; }
  double  getPermZ(int& i)  { return zPerm[i]; }
  int getIndex(int i, int j, int k )    { return i + j * nx     + k * nx     * ny; }
  int getPermIndex(int i, int j, int k) { return i + j * nxPerm 
	                                           + k * nxPerm * nyPerm;}
  int toIPerm(int i) { return iPerm[i]; }
  int toJPerm(int j) { return jPerm[j]; }
  int toKPerm(int k) { return kPerm[k]; } 
	  
  GridType getGridType() { return gridType; }

  // --- write nx, ny, nz, x,y,z to file, for later transport part ---
  void writeXYZ( Unit unit );
  void PrintGrid( ostream & );

 private:
  bool debug;

  Unit unit;
  GridType gridType;
  char   *directory;

  // pressure grid
  int nx, ny, nz, nNode;
  double *x, *y, *z, dim[3];    // size: nx, ny, nz;
  double *dx, *dy, *dz;         // size: nx-1, ny-1, nz-1;
  double *bdx, *bdy, *bdz;      // size: nx, ny, nz;

  // perm grid
  int    nxPerm, nyPerm, nzPerm, nNodePerm;
  int    *iPerm, *jPerm, *kPerm;
  double *xPerm, *yPerm, *zPerm;

  // --- For grid scope data
  void readData( ifstream &is );
  void initialize();

  // --- Calculate geometry data (x,dx,bdx,dim[0].)
  void setupGrid();
  void setupPermGrid();
  
  // from Field unit to Metric unit
  void convertUnit();
};
#endif
