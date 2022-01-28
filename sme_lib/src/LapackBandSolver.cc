/*
 * File: LapackBandSolver.cc
 * ----------------------------------
 * Implementation for LapackBandSolver 
 */
#include <math.h>
#include <time.h>
#include "Eqn.h" 
#include "LapackBandSolver.h" 
 
/* --- Public Methods ----------------------------------- */
// --- Constructors and destructor --------------
LapackBandSolver::LapackBandSolver(const Eqn &e) : Solver(e)
{ 
  if(nx != 1 && ny != 1) { 
    if(nz == 1) {             // xy (always use nx, and no bandTransform)
      bandSize = nx;   
    } else {                  // 3D (use ny*nz with bandTransform, or nx*ny with no) 
      bandSize = (nx > nz) ? ny*nz : nx*ny;
    }
  } else {                    // 1D, 2D(xz)
    int nxyzMax = nx;  
    if(ny>nxyzMax)  nxyzMax = ny;
    if(nz>nxyzMax)  nxyzMax = nz;    // the max(nx, ny, nz)
    bandSize = nx*ny*nz/nxyzMax;     // the product of two smallest number
  }

  bandSwitch = (nz > 1 && nx > nz);  // band Tranform criterian
 
  aHeight = 3 * bandSize + 1;
  A = new double[ aHeight * nA ];
  ipiv = new int[ nA ];
}

LapackBandSolver::~LapackBandSolver()
{
  delete[] A;   delete[] ipiv;  
}

// --- solve AX=B --------------------------------
void LapackBandSolver::solve(bool doLU, const Eqn &e, int nRHS, double *RHS)
{
  int i, j, k, ir;
  double startT, endT, *RHS1;  // RHS1: tmp storage for one RHS
 
  if( doLU ) { 
    // --- A is just assigned or recalculated, need do LU Decomp to A - 
    //cout << "NOW, DO LU DECOMP TO MATRIX A --------------- " << endl;
    // --- first setup A from Eqn --------------
    if( bandSwitch ) setupA_Transform(e);      // band transform
    else setupA(e);                            // no band transform

    // --- do LU factorization ------------
    startT = clock();
    int succ = DGBTRF_C(nA, nA, bandSize, bandSize, A, aHeight, ipiv);
    if( succ != 0) {
        cerr << "ERROR in LU ----------------\n";
        cerr << "INFO = " << succ << endl;
        cerr << "Matrix dimension = " << nA << " x " << nA << endl;
    }
    endT = clock();
    //cerr << "LUDecomp time: "<< (endT-startT)/CLOCKS_PER_SEC << " sec" << endl;
  }

  // --- if bandTransform is performed, 
  // --- for each RHS, change index order from ijk to kji ---
  if( bandSwitch ) {   
      RHS1 = new double[nA];   
      for(ir=0; ir<nRHS; ir++) {  // transform one RHS
      for(k=0; k<nz; k++) {
      for(j=0; j<ny; j++) {
      for(i=0; i<nx; i++) { 
          RHS1[k+(j+i*ny)*nz] = RHS[(i+(j+k*ny)*nx)+ir*nA];   // ijk to kji
      }
      }
      }
          for(i=0; i<nA; i++) RHS[i+ir*nA] = RHS1[i]; // put it back in place
      }
  }

  // --- back solve -----------
  DGBTRS_C('N', nA, bandSize, bandSize, nRHS, A, aHeight, ipiv, RHS, nA);

  // --- if bandTransform is performed, 
  // --- for each X, change index order from kji to ijk ---
  if( bandSwitch ) {
      for(ir = 0; ir < nRHS; ir++) {  // for each X
      for(k = 0; k < nz; k++) {
      for(j = 0; j < ny; j++) {
      for(i = 0; i < nx; i++) { 
          RHS1[i+(j+k*ny)*nx] = RHS[k+(j+i*ny)*nz+ir*nA];  // from kji to ijk
      }
      }
      }
          for(i=0; i<nA; i++) RHS[i+ir*nA] = RHS1[i];   // put it back in place
      }
      delete[] RHS1;  // remove tmp storage
  }
}

/* --- Protected Methods -------------------------------------- */
// --- setup A without band transform -------
// --- refer to Lapack menu for "dgbtrf" for A's storage format ---
void LapackBandSolver::setupA(const Eqn &e)
{
  int i;
  // --- Note: A is column wise in Fortran!!!!!!!  --- 

  for(i = 0; i < nA * aHeight; i++) A[i] = 0.0;        // all zero
  
  for(i = 0; i < nA; i++) 
      A[2 * bandSize + i * aHeight] = e.diagA[i];      // assign diagonal

  if( nx > 1) {                                 // assign off-D X+ and X-
      for(i = 1; i < nA    ; i++) 
          A[2 * bandSize - 1 + i * aHeight] = e.offDAXp[i-1];
      for(i = 0; i < nA - 1; i++) 
          A[2 * bandSize + 1 + i * aHeight] = e.offDAXn[i  ];
  }
  
  if( ny > 1) {                                   // assign off-D Y+ and Y-
      for(i = nx; i < nA     ; i++) 
          A[2 * bandSize -nx + i * aHeight] = e.offDAYp[i-nx];
      for(i =  0; i < nA - nx; i++) 
          A[2 * bandSize + nx + i * aHeight] = e.offDAYn[i];
  }
  
  if( nz > 1) {                                   // assign off-D Z+ and Z-
      int nxy = nx*ny;
      for(i = nxy; i < nA; i++)       
          A[2 * bandSize - nxy + i * aHeight] = e.offDAZp[i-nxy];
      for(i = 0;   i < nA - nxy; i++) 
          A[2 * bandSize + nxy + i * aHeight] = e.offDAZn[i];
  }
}

// --- setup A with band transform ------------------ 
// --- now, I bands are at ny*nz, J bands are at nz, K bands are at 1 ---
// --- refer to Lapack menu for "dgbtrf" for A's storage format ---
void LapackBandSolver::setupA_Transform(const Eqn &e)
{
  int i, j, k, nxy=nx*ny, nyz=ny*nz;

  // --- Note: A is column wise in Fortran!!!!!!!  --- 
  for(i=0; i<nA*aHeight; i++) A[i]=0.0;         // all zero

  for(k=0; k<nz; k++) {                         // assign diagonal 
  for(j=0; j<ny; j++) {                         // from ijk to kji
  for(i=0; i<nx; i++) {
    A[2*bandSize+(k+(j+i*ny)*nz)*aHeight]=e.diagA[(i+(j+k*ny)*nx)];
  }}}

  if(nx>1) {                                    // assign off-D X+ and X-
    for(k=0; k<nz; k++) {                       // from ijk to kji
    for(j=0; j<ny; j++) {                       // put I bands at ny*nz+- 
      for(i=1; i<nx; i++) 
        A[2*bandSize-nyz+(k+(j+i*ny)*nz)*aHeight]=e.offDAXp[(i+(j+k*ny)*nx)-1];
      for(i=0; i<nx-1; i++) 
        A[2*bandSize+nyz+(k+(j+i*ny)*nz)*aHeight]=e.offDAXn[(i+(j+k*ny)*nx)];
    }}
  }
  if(ny>1) {                                    // assign off-D Y+ and Y-
    for(k=0; k<nz; k++) {                       // from ijk to kji
    for(i=0; i<nx; i++) {                       // put J bands at nz+-
      for(j=1; j<ny; j++) 
        A[2*bandSize-nz+(k+(j+i*ny)*nz)*aHeight]=e.offDAYp[(i+(j+k*ny)*nx)-nx];
      for(j=0; j<ny-1; j++) 
        A[2*bandSize+nz+(k+(j+i*ny)*nz)*aHeight]=e.offDAYn[(i+(j+k*ny)*nx)];
    }}
  }
  if(nz>1) {                                    // assign off-D Z+ and Z-
    for(j=0; j<ny; j++) {                       // from ijk to kji
    for(i=0; i<nx; i++) {                       // put K bands at 1+-
      for(k=1; k<nz; k++) 
        A[2*bandSize-1+(k+(j+i*ny)*nz)*aHeight]=e.offDAZp[(i+(j+k*ny)*nx)-nxy];
      for(k=0; k<nz-1; k++) 
        A[2*bandSize+1+(k+(j+i*ny)*nz)*aHeight]=e.offDAZn[(i+(j+k*ny)*nx)];
    }}
  }
}

/* --- Private Methods ----------------------------- */
// --- C++ requirement to call Fortran subroutine ----
extern "C"
{
  void dgbtrf_(int *m, int *n, int *kl, int *ku, double *a, int *lda,
                 int *ipiv, int *info);
  void dgbtrs_(char *tran, int *n, int *kl, int *ku, int *nRHS, 
               double *a, int *lda, int *ipiv, 
               double *b, int *ldb, int *info);
};

// --- Wrappers to Lapack Fortran subroutines ------
// --- band matrix LU decomposition --------------
int LapackBandSolver::DGBTRF_C(int m, int n, int kl, int ku, double *a, 
                              int lda, int *ipiv)
{
  int info = -1;
  dgbtrf_(&m, &n, &kl, &ku, a, &lda, ipiv, &info);
  return info;
}

// --- band matrix back solve using band LU results -------------
int LapackBandSolver::DGBTRS_C(char tran, int n, int kl, int ku, int nRHS, 
                               double *a, int lda, int *ipiv, 
                               double *b, int ldb)
{
  int info = -1;
  dgbtrs_(&tran, &n, &kl, &ku, &nRHS, a, &lda, ipiv, b, &ldb, &info);
  return info;
}

