/*
 * File: LapackSymBandSolver.cc
 * ----------------------------------
 * Implementation for LapackSymBandSolver 
 */
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "Eqn.h"
#include "LapackSymBandSolver.h"

/* --- Public Methods -------------------------------- */
/* --- Constructors and destructor --- */
LapackSymBandSolver::LapackSymBandSolver(const Eqn &e) : Solver(e)
{ 
  if(nx!=1&&ny!=1) { 
    if(nz==1) {        // xy (always use nx, and no bandTransform)
      bandSize = nx;   
    } else {           // 3D (use ny*nz with bandTransform, or nx*ny with no) 
      bandSize = (nx>nz)? ny*nz : nx*ny;
    }
  } else {                    // 1D, 2D(xz)
    int nxyzMax = nx;  
    if(ny>nxyzMax)  nxyzMax = ny;
    if(nz>nxyzMax)  nxyzMax = nz;    // the max(nx, ny, nz)
    bandSize = nx*ny*nz/nxyzMax;     // the product of two smallest number
  }
  
  bandSwitch = (nz>1&&nx>nz);        // band Tranform criterian
  
  aHeight = bandSize+1;
  A = new double[aHeight*nA];
}

LapackSymBandSolver::~LapackSymBandSolver()
{
  delete[] A;  
}

// --- solve AX=B ---------------------------
void LapackSymBandSolver::solve(bool doLU, const Eqn &e, int nRHS, double *RHS)
{
  int i, j, k, ir; 
  double startT, endT, *RHS1;   // RHS1: tmp storage for one RHS
 
  if(doLU) { 
    // --- A is just assigned or recalculated, need do LU Decomp to A - 
    cout << "NOW, DO LU DECOMP TO MATRIX A ---------------- " << endl;
    
    // --- symmetric condition checking ---
    cout << "CHECK SYMMETRIC CONDITION OF MATRIX A -------- " << endl;
    if(nx>1) {
      for(i=0; i<nA-1; i++) if(fabs(e.offDAXp[i]-e.offDAXn[i])>1e-10) { 
	cerr << "ERROR,non-symmetric in X+-: " << i << endl;  exit(0); }
    }
    if(ny>1) { 
      for(i=0; i<nA-nx; i++) if(fabs(e.offDAYp[i]-e.offDAYn[i])>1e-10) { 
	cerr << "ERROR, non-symmetric in Y+-: " << i << endl; exit(0); }
    }
    if(nz>1) {
      int nxy = nx*ny;
      for(i=0; i<nA-nxy; i++) if(fabs(e.offDAZp[i]-e.offDAZn[i])>1e-10) {
	cerr << "ERROR, non-symmetric in Z+-: " << i << endl; exit(0); }
    }

    // --- first setup A from Eqn --------------
    if(bandSwitch) setupA_Transform(e);      // band transform
    else setupA(e);                          // no band transform

    // --- do LU factorization ------------
    startT = clock();
    int succ = DPBTRF_C('U', nA, bandSize, A, aHeight);
    if(succ!=0) {
      cerr << "ERROR in LU (possible not positive defined)---"<< succ <<"-\n";
      exit(0);
    }
    endT = clock();
    cerr << "LUDecomp time: "<< (endT-startT)/CLOCKS_PER_SEC << " sec" << endl;
  }

  // --- add a "-" to RHS ro compensate adding a "-" to A --- 
  for(i=0; i<nA*nRHS; i++) RHS[i] = -RHS[i];

  // --- if bandTransform is performed, 
  // --- for each RHS, change index order from ijk to kji ---
  if(bandSwitch) {
    RHS1 = new double[nA];   
    for(ir=0; ir<nRHS; ir++) {  // transform one RHS
      for(k=0; k<nz; k++) {
      for(j=0; j<ny; j++) {
      for(i=0; i<nx; i++) { 
	RHS1[k+(j+i*ny)*nz] = RHS[(i+(j+k*ny)*nx)+ir*nA];   // ijk to kji
      }}}
      for(i=0; i<nA; i++) RHS[i+ir*nA] = RHS1[i]; // put it back in place
    }
  }

  // --- back solve -----------
  DPBTRS_C('U', nA, bandSize, nRHS, A, aHeight, RHS, nA);

  // --- if bandTransform is performed, 
  // --- for each X, change index order from kji to ijk ---
  if(bandSwitch) {
    for(ir=0; ir<nRHS; ir++) {  // for each X 
      for(k=0; k<nz; k++) {
      for(j=0; j<ny; j++) {
      for(i=0; i<nx; i++) { 
	RHS1[i+(j+k*ny)*nx] = RHS[k+(j+i*ny)*nz+ir*nA];  // from kji to ijk
      }}}
      for(i=0; i<nA; i++) RHS[i+ir*nA] = RHS1[i];   // put it back in place
    }
  }
   
  if(bandSwitch) delete[] RHS1; // remove tmp storage
}

/* --- Protected Methods -------------------------------------- */
// --- setup A without band transform -------
// --- refer to Lapack menu for "dpbtrf" for A's storage format ---
// --- due to symmetric contidion only upper half is assigned ---
// --- A is added a "-" to satisfy the positive define condition ---
void LapackSymBandSolver::setupA(const Eqn &e)
{
  int i;
  // --- Note: A is column wise in Fortran!!!!!!!  --- 
  for(i=0; i<nA*aHeight; i++) A[i]=0.0;           // all zero

  for(i=0; i<nA; i++) A[bandSize+i*aHeight]=-e.diagA[i]; // assign diagonal 
    
  if(nx>1) {                                      // assign off-D X+
    for(i=1; i<nA; i++) A[bandSize-1+i*aHeight]=-e.offDAXp[i-1];
  }
  if(ny>1) {                                      // assign off-D Y+
    for(i=nx; i<nA; i++) A[bandSize-nx+i*aHeight]=-e.offDAYp[i-nx];
  }
  if(nz>1) {                                      // assign off-D Z+
    int nxy = nx*ny;
    for(i=nxy; i<nA; i++) A[bandSize-nxy+i*aHeight]=-e.offDAZp[i-nxy];
  }
}

// --- setup A with band transform -------
// --- now, I bands are at ny*nz, J bands are at nz, K bands are at 1 ---
// --- refer to Lapack menu for "dpbtrf" for A's storage format ---
// --- due to symmetric contidion only upper half is assigned ---
// --- A is added a "-" to satisfy the positive define condition ---
void LapackSymBandSolver::setupA_Transform(const Eqn &e)
{
  int i, j, k, nxy=nx*ny, nyz=ny*nz;

  // --- Note: A is column wise in Fortran!!!!!!!  --- 
  for(i=0; i<nA*aHeight; i++) A[i]=0.0;         // all zero

  for(k=0; k<nz; k++) {                         // assign diagonal 
  for(j=0; j<ny; j++) {                         // from ijk to kji
  for(i=0; i<nx; i++) {
    A[2*bandSize+(k+(j+i*ny)*nz)*aHeight]=e.diagA[(i+(j+k*ny)*nx)];
  }}}
    
  if(nx>1) {                                    // assign off-D X+
    for(k=0; k<nz; k++) {                       // from ijk to kji
    for(j=0; j<ny; j++) {                       // put I+ band at ny*nz+ 
    for(i=1; i<nx; i++) {
      A[bandSize-nyz+(k+(j+i*ny)*nz)*aHeight]=-e.offDAXp[(i+(j+k*ny)*nx)-1];
    }}}
  }
  if(ny>1) {                                    // assign off-D Y+ 
    for(k=0; k<nz; k++) {                       // from ijk to kji
    for(i=0; i<nx; i++) {                       // put J+ band at nz+
    for(j=1; j<ny; j++) {
      A[bandSize-nz+(k+(j+i*ny)*nz)*aHeight]=-e.offDAYp[(i+(j+k*ny)*nx)-nx];
    }}}
  }
  if(nz>1) {                                    // assign off-D Z+
    for(j=0; j<ny; j++) {                       // from ijk to kji
    for(i=0; i<nx; i++) {                       // put K+ band at 1+
    for(k=1; k<nz; k++) { 
      A[bandSize-1+(k+(j+i*ny)*nz)*aHeight]=-e.offDAZp[(i+(j+k*ny)*nx)-nxy];
    }}}
  }
}

/* --- Private Methods ----------------------------- */
// --- C++ requirement to call Fortran subroutine ----
extern "C"
{
  void dpbtrf_(char *uplo, int *n, int *kd, double *a, int *lda, 
		       int *info);
  void dpbtrs_(char *uplo, int *n, int *kd, int *nRHS, double *a, 
		       int *lda, double *b, int *ldb, int *info);
};

// --- Wrappers to Lapack Fortran subroutines ------
// --- symmetric band matrix Cholesky decomposition ---------
int LapackSymBandSolver::DPBTRF_C(char uplo, int n, int kd, double *a, 
				  int lda)
{
  int info = -1;
  dpbtrf_(&uplo, &n, &kd, a, &lda, &info);
  return info;
}

// --- symmetric band matrix back solve using sym band LU results --------
int LapackSymBandSolver::DPBTRS_C(char uplo, int n, int kd, int nRHS, 
			       double *a, int lda, double *b, int ldb)
{
  int info = -1;
  dpbtrs_(&uplo, &n, &kd, &nRHS, a, &lda, b, &ldb, &info);
  return info;
}

