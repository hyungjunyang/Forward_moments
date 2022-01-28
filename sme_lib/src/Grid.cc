/*
 * File: Grid.cc
 * ----------------------------------
 * Implementation for Grid class
 */
                              
#include "Grid.h"           
                               
//===== ( Grid( ifstream &is, Unit u, char *Dir ) ) =====
Grid::Grid( ifstream &is, Unit u, char *Dir ) {
  debug = true;

  if(debug) cout << "Grid::Grid()\n";
  
  directory = Dir;
  unit = u;
  readData(is);  // all initialization are inside
  
  setupGrid();
  setupPermGrid();
  if( unit == FIELD) {
      convertUnit();  // from field to metric
      //for(int i = 0; i < nx; i++) {
      //    cout << x[i]<< "  " << bdx[i] << endl; 
      //}  
      //exit(0);      
  }
  writeXYZ( unit );
}

//===== ( ~Grid() ) =====
  Grid::~Grid() {
  if(debug) cout << "Grid::~Grid()" << endl;
  delete[] x;    delete[] y;    delete[] z;
  delete[] dx;   delete[] dy;   delete[] dz;
  delete[] bdx;  delete[] bdy;   delete[] bdz;
  delete[] xPerm; delete[] yPerm; delete[] zPerm;
  delete[] iPerm; delete[] jPerm; delete[] kPerm;
}

//===== ( readData( ifstream &is  ) =====
  void Grid::readData( ifstream &is ) {
  int i;
  double tmp; 

  // --- grid geometry data ---   
  Junk(is); is >> i; gridType = (GridType) i;
  Junk(is); is >> nx >> ny >> nz; 
  if(debug) {
     cout << i << endl;
     cout << nx << ' ' 
          << ny << ' ' 
          << nz << endl;
     //exit(0);
  }
  
  nNode = nx * ny * nz;
  
  initialize();    // initialize grid scope data
  if(debug) {
    cout << " finish initialization" << endl; 
  }
  // point distributed  (read dx,  dy,  dz)  
  if( gridType == POINT_DIS) {
      tmp = flagReading(is, nx-1, dx); if(nx == 1) bdx[0] = tmp;
      tmp = flagReading(is, ny-1, dy); if(ny == 1) bdy[0] = tmp;
      tmp = flagReading(is, nz-1, dz); if(nz == 1) bdz[0] = tmp;
      // Note: nx = 1, no dx, but still has one bdx (A special Condition)
      // the input at this condition is the physical dimension or bdx[0]
  }
  // cell centered      (read bdx, bdy, bdz)
  else if( gridType == CELL_CEN ) {
      flagReading(is, nx, bdx);
      flagReading(is, ny, bdy);
      flagReading(is, nz, bdz);
  }
  else {
      cerr << "Error in Grid Type \n";
  }
}

//===== ( initialize() ) =======================================================
void Grid::initialize() {
  x = new double[nx];
  y = new double[ny];
  z = new double[nz];
  
  bdx = new double[nx];
  bdy = new double[ny];
  bdz = new double[nz];

  dx  = (nx > 1) ? new double[nx-1] : NULL;  // NULL is OK for delete..
  dy  = (ny > 1) ? new double[ny-1] : NULL;
  dz  = (nz > 1) ? new double[nz-1] : NULL;
}

//===== ( setupGrid() ) ========================================================
void Grid::setupGrid() {
  int i, j, k;
  // --- calculate bdx.. or dx... --------------
  if( gridType == POINT_DIS ) {  // point distributed grid
    if( nx > 1) { bdx[0] = dx[0]/2.0; bdx[nx-1] = dx[nx-2]/2.0; } // boundary cells
    if( ny > 1) { bdy[0] = dy[0]/2.0; bdy[ny-1] = dy[ny-2]/2.0; }
    if( nz > 1) { bdz[0] = dz[0]/2.0; bdz[nz-1] = dz[nz-2]/2.0; }
    for( i = 1; i < nx - 1; i++) bdx[i] = (dx[i-1] + dx[i])/2.0;  // internal cells
    for( j = 1; j < ny - 1; j++) bdy[j] = (dy[j-1] + dy[j])/2.0;
    for( k = 1; k < nz - 1; k++) bdz[k] = (dz[k-1] + dz[k])/2.0;
  } 
  else {                   // cell centered grid
    for(i = 0; i < nx - 1; i++) dx[i] = (bdx[i] + bdx[i+1])/2.0;
    for(j = 0; j < ny - 1; j++) dy[j] = (bdy[j] + bdy[j+1])/2.0;
    for(k = 0; k < nz - 1; k++) dz[k] = (bdz[k] + bdz[k+1])/2.0;
  }

  // --- calculate x, y, z ------------------------
  if(gridType == POINT_DIS) { 
     x[0] = y[0] = z[0] = 0.0 ; 
  }
  else { 
     x[0] = bdx[0]/2.0; y[0]=bdy[0]/2.0; z[0]=bdz[0]/2.0; 
  }

  for(i = 0; i < nx - 1; i++) x[i+1] = x[i] + dx[i];
  for(j = 0; j < ny - 1; j++) y[j+1] = y[j] + dy[j];
  for(k = 0; k < nz - 1; k++) z[k+1] = z[k] + dz[k];

  // --- calculate physical dimension -----
  dim[0] = dim[1] = dim[2] = 0.0;
  for(i = 0; i < nx; i++) dim[0] += bdx[i];
  for(j = 0; j < ny; j++) dim[1] += bdy[j];
  for(k = 0; k < nz; k++) dim[2] += bdz[k];
}

//===== ( setupPermGrid() ) ====================================================
void Grid::setupPermGrid(){
  if( gridType == POINT_DIS) {
      nxPerm = nx * 2 - 1;
      nyPerm = ny * 2 - 1;
      nzPerm = nz * 2 - 1;
  }
  else if( gridType == CELL_CEN ) {
      if( nx > 1 ) nxPerm = nx * 2 + 1;
          else     nxPerm = nx;
      if( ny > 1 ) nyPerm = ny * 2 + 1;
          else     nyPerm = ny;
      if( nz > 1 ) nzPerm = nz * 2 + 1;
          else     nzPerm = nz;
  } 
  else {
      cerr << "Error in Grid Type \n";
  }
  nNodePerm = nxPerm * nyPerm * nzPerm;
  xPerm     = new double[ nxPerm ];
  yPerm     = new double[ nyPerm ];
  zPerm     = new double[ nzPerm ];
  iPerm     = new int   [ nx ];
  jPerm     = new int   [ ny ];
  kPerm     = new int   [ nz ];

  int i, j, k;
  if( gridType == POINT_DIS ) {
      for( i = 0; i < nx; i++ ) {
	          iPerm[i] = 2 * i;
           xPerm[ iPerm[i] ] = x[i];
      }
      for( i = 0; i < ny; i++ ) {
	          jPerm[i] = 2 * i;
           yPerm[ jPerm[i] ] = y[i];
      }
      for( i = 0; i < nz; i++ ) {
	          kPerm[i] = 2 * i;
           zPerm[ kPerm[i] ] = z[i];
      }
      
      if( nx > 1 ) for( i = 1; i < 2 * nx -1; i = i + 2) {
           xPerm[i] = ( xPerm[i-1] + xPerm[i + 1] )/2. ;
      }
      if( ny > 1 ) for( i = 1; i < 2 * ny -1; i = i + 2) {
           yPerm[i] = ( yPerm[i-1] + yPerm[i + 1] )/2. ;
      }
      if( nz > 1 ) for( i = 1; i < 2 * nz -1; i = i + 2) {
           zPerm[i] = ( zPerm[i-1] + zPerm[i + 1] )/2. ;
      }
  }
  else if( gridType == CELL_CEN ) {
      if (nx > 1) {
	 xPerm[0] = 0.;
         for( i = 0; i < nx; i++ ) {
	   iPerm[i] = 2 * i + 1;
           xPerm[ iPerm[i] ]     = x[i];
	   xPerm[ iPerm[i] + 1 ] = xPerm[ iPerm[i] - 1 ] + bdx[i];
	 }
      }
      else {
	   iPerm[0] = 0; 
           xPerm[ iPerm[0] ] = x[0];
      }
      if (ny > 1) {
	 yPerm[0] = 0; 
         for( i = 0; i < ny; i++ ) {
	   jPerm[i] = 2 * i + 1;
           yPerm[ jPerm[i] ]     = y[i];
	   yPerm[ jPerm[i] + 1 ] = yPerm[ jPerm[i] - 1 ] + bdy[i];
         }
      }
      else {
	   jPerm[0] = 0; 
           yPerm[ jPerm[0] ] = y[0];
      }
      if (nz > 1) {
         zPerm[0] = 0.;
         for( i = 0; i < nz; i++ ) {
	   kPerm[i] = 2 * i + 1;
           zPerm[ kPerm[i] ]     = z[i];
	   zPerm[ kPerm[i] + 1 ] = zPerm[ kPerm[i] - 1 ] + bdz[i];
         }
      }
      else {
	   kPerm[0] = 0; 
           zPerm[ kPerm[0] ] = z[0];
      }
  }
  else {
      cerr << "Error in Grid Type \n"; 
  }
}

//===== ( PrintGrid( ostream &os ) ) =====
void Grid::PrintGrid( ostream &os ) {
  os << endl;
  os << "===== Grid Geometry =====" << endl;
  (gridType == POINT_DIS)? os << "Point Distributed Grid" << endl 
    :    os << "Cell Centered Grid" << endl;
  os << "Grid Dimensions: (" << nx << ", " 
                             << ny << ", " 
                             << nz << ")\n";
  os << "Physical dimensions : (" << dim[0] << ", " 
                                  << dim[1] << ", " 
                                  << dim[2] << ")\n";
}

//===== ( writeXYZ(Unit unit) ) ================================================
void Grid::writeXYZ(Unit unit) {
  int i;
  char outfile[20],outgrid[20];
  sprintf( outfile, "%sxyz.out", directory );
  sprintf( outgrid, "%sgrid.out", directory );
  ofstream os(outfile, ios::out);
  ofstream os1(outgrid, ios::out);

  for(int j = 0; j < nyPerm; j = j + 2) {
      for( i = 0; i < nxPerm; i = i + 2) {
           os1 << xPerm[i]<< "  " << yPerm[j] << endl; 
      }            
      os1 << endl;
  }
  for(int j = 0; j < nyPerm; j = j + 2) {
      for( i = 0; i < nxPerm; i = i + 2) {
           os1 << yPerm[j]<< "  " << xPerm[i] << endl; 
      }            
      os1 << endl;
  }
  os1.close();
  
  os1.close();
  
  os << nx << "  " << ny << "  " << nz << endl;
  
  // --- write dim[0], dim[1], dim[2] ---
  (unit==FIELD)? 
    os << dim[0]/ft_m << "  " << dim[1]/ft_m << "  " << dim[2]/ft_m : 
    os << dim[0]      << "  " << dim[1] << "  " << dim[2];
  os << endl << endl;

  for(i=0; i<nx; i++) {
    (unit==FIELD)? os << x[i]/ft_m << "  " : os << x[i] << "  "; 
    if((i+1)%6==0) os << endl;
  }
  os << endl << endl;
  
  for(i=0; i<ny; i++) {
    (unit==FIELD)? os << y[i]/ft_m << "  " : os << y[i] << "  "; 
    if((i+1)%6==0) os << endl;
  }
  os << endl << endl;
  
  for(i=0; i<nz; i++) {
    (unit==FIELD)? os << z[i]/ft_m << "  " : os << z[i] << "  "; 
    if((i+1)%6==0) os << endl;
  }
  os << endl << endl;

  
  // --- write nx, ny, nz ---
  os << "nx = " << nx << "  " 
     << "ny = " << ny << "  " 
     << "nz = " << nz << endl;
  
  // --- write dim[0], dim[1], dim[2] ---
  os << "Lx "  << "Ly "  << "Lz = : ";
  (unit == FIELD)? 
    os << dim[0]/ft_m << "  " << dim[1]/ft_m << "  " << dim[2]/ft_m : 
    os << dim[0] << "  " << dim[1] << "  " << dim[2];
  os << endl << endl;

  // --- write x[], y[], z[] ---------
  os << "x[i] = 0 ... " << nx-1<< endl; 
  for( i = 0; i < nx; i++) {
       (unit == FIELD)? os << x[i]/ft_m << "  " : os << x[i] << "  "; 
       if( (i+1) % 9 == 0) os << endl;
  }            
  os << endl << endl;
  
  os << "xPerm[i] = 0 ... " << nxPerm-1<< endl; 
  for( i = 0; i < nxPerm; i++) {
       (unit == FIELD)? os << xPerm[i]/ft_m << "  " : os << xPerm[i] << "  "; 
       if( (i+1) %9 == 0) os << endl;
  }            
  os << endl << endl;

  os << "y[j] = 0 ... " << ny-1<< endl; 
  for( i = 0; i < ny; i++) {
       (unit == FIELD)? os << y[i]/ft_m << "  " : os << y[i] << "  "; 
       if( (i+1) %9 == 0) os << endl;
  }
  os << endl << endl;
  os << "yPerm[j] = 0 ... " << nyPerm-1<< endl; 
  for( i = 0; i < nyPerm; i++) {
       (unit == FIELD)? os << yPerm[i]/ft_m << "  " : os << yPerm[i] << "  "; 
       if( (i+1) %9 == 0) os << endl;
  }
  os << endl << endl;
  
  os << "z[k] = 0 ... " << nz - 1<< endl; 
  for( i = 0; i < nz; i++) {
       (unit == FIELD)? os << z[i]/ft_m << "  " : os << z[i] << "  "; 
       if( (i+1) %9 ==0) os << endl;
  }
  os << endl << endl;
  os << "zPerm[k] = 0 ... " << nzPerm - 1<< endl; 
  for( i = 0; i < nzPerm; i++) {
       (unit == FIELD)? os << zPerm[i]/ft_m << "  " : os << zPerm[i] << "  "; 
       if( (i+1) %9 ==0) os << endl;
  }
  os << endl << endl;

  // --- write dx[], dy[], dz[] ---------
  os << "dx[i] = 0 ... " << nx - 2<< endl; 
  for( i = 0; i < nx - 1; i++) {
       (unit == FIELD)? os << dx[i]/ft_m << "  " : os << dx[i] << "  "; 
       if( (i+1) %6 == 0) os << endl;
  }            
  os << endl << endl;
  os << "dy[j] = 0 ... " << ny - 2<< endl; 
  for( i = 0; i < ny - 1; i++) {
       (unit == FIELD)? os << dy[i]/ft_m << "  " : os << dy[i] << "  "; 
       if( (i+1) %6 == 0) os << endl;
  }
  os << endl << endl;
  os << "dz[k] = 0 ... " << nz - 2<< endl; 
  for( i = 0; i < nz - 1; i++) {
       (unit == FIELD)? os << dz[i]/ft_m << "  " : os << dz[i] << "  "; 
       if( (i+1) %6 ==0) os << endl;
  }
  os << endl << endl;
  
  // --- write bdx[], bdy[], bdz[] ---------
  os << "bdx[i] = 0 ... " << nx - 1<< endl; 
  for( i = 0; i < nx ; i++) {
       (unit == FIELD)? os << bdx[i]/ft_m << "  " : os << bdx[i] << "  "; 
       if( (i+1) %9 == 0) os << endl;
  }
  os << endl << endl;
  os << "bdy[j] = 0 ... " << ny - 1<< endl; 
  for( i = 0; i < ny ; i++) {
       (unit == FIELD)? os << bdy[i]/ft_m << "  " : os << bdy[i] << "  "; 
       if( (i+1) %9 == 0) os << endl;
  }
  os << endl << endl;
  os << "bdz[k] = 0 ... " << nz - 1<< endl; 
  for( i = 0; i < nz ; i++) {
       (unit == FIELD)? os << bdz[i]/ft_m << "  " : os << bdz[i] << "  "; 
       if( (i+1) %9 ==0) os << endl;
  }
  os << endl << endl;
  
  os <<"Xperm used in code" << endl;
  for( i = 0; i < nxPerm; i++) {
       os << "i = "    << i << ' ' 
	  << "xperm = "<< xPerm[i] << "  "; 
       if( (i+1) %2 ==0) os << endl;

  }            
  os << endl << endl;
 
  os.close();
}

//===== ( convertUnit() ) ======================================================
void Grid::convertUnit() {
  int i, j, k;
  for(i = 0; i < nx; i++) { x[i] *= ft_m; bdx[i] *= ft_m; }
  for(j = 0; j < ny; j++) { y[j] *= ft_m; bdy[j] *= ft_m; }
  for(k = 0; k < nz; k++) { z[k] *= ft_m; bdz[k] *= ft_m; }

  for(i = 0; i < nx - 1; i++) dx[i]  *= ft_m; 
  for(j = 0; j < ny - 1; j++) dy[j]  *= ft_m; 
  for(k = 0; k < nz - 1; k++) dz[k]  *= ft_m; 
  for(i = 0; i < 3     ; i++) dim[i] *= ft_m;

  for(i = 0; i < nxPerm; i++) { xPerm[i] *= ft_m; }
  for(j = 0; j < nyPerm; j++) { yPerm[j] *= ft_m; }
  for(k = 0; k < nzPerm; k++) { zPerm[k] *= ft_m; }
}
