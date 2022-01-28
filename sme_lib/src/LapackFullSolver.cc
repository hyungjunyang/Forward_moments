/*
 * File: LapackFullSolver.cc
 * ----------------------------------
 * Implementation for LapackFullSolver class
 */
#include "LapackFullSolver.h"
#include <math.h>
#include <time.h>
#include "Eqn.h"

/* --- Public Methods ----------------------------------------- */
// --- Constructors and destructor ---------------------------
LapackFullSolver::LapackFullSolver(const Eqn &e) : Solver(e)
{ 
  A = new double[nA*nA];
  ipiv = new int[nA];
}

LapackFullSolver::~LapackFullSolver()
{
  delete[] A;  delete[] ipiv;
}

// --- solve AX=B ----------------------------------
void LapackFullSolver::solve(bool doLU, const Eqn &e, int nRHS, double *RHS)
{
  double startT, endT;
 
  if(doLU) { 
    // --- A is just assigned or recalculated, need do LU Decomp to A - 
    cout << "NOW, DO LU DECOMP TO MATRIX A ---------------- " << endl;
    setupA(e);      // first setup A from Eqn
    
    // --- do LU factorization ------------
    startT = clock();
    int succ = DGETRF_C(nA, nA, A, nA, ipiv);
    if(succ!=0) cerr << "ERROR in LU ----------------\n";
    endT = clock();
    cerr << "LUDecomp time: "<<(endT-startT)/CLOCKS_PER_SEC << " sec" << endl;
  }

  // --- back solve -----------
  DGETRS_C('N', nA, nRHS, A, nA, ipiv, RHS, nA);
}

/* --- Protected Methods -------------------------------- */
// --- setup matrix A -----------
void LapackFullSolver::setupA(const Eqn &e)
{
  int i;
  // --- Note: A is column wise in Fortran!!!!!!!  --- 
  for(i=0; i<nA*nA; i++) A[i] = 0.0;           // all zero
  
  for(i=0; i<nA; i++) A[i+i*nA]=e.diagA[i];    // assign diagonal 
  if(nx>1) {                                   // assign off-D X+ and X-
    for(i=1; i<nA; i++) A[i-1+i*nA]=e.offDAXp[i-1];
    for(i=0; i<nA-1; i++) A[i+1+i*nA]=e.offDAXn[i];
  }
  if(ny>1) {                                   // assign off-D Y+ and Y-
    for(i=nx; i<nA; i++) A[i-nx+i*nA]=e.offDAYp[i-nx];
    for(i=0; i<nA-nx; i++) A[i+nx+i*nA]=e.offDAYn[i];
  }
  if(nz>1) {                                   // assign off-D Z+ and Z-
    int nxy = nx*ny;
    for(i=nxy; i<nA; i++) A[i-nxy+i*nA]=e.offDAZp[i-nxy];
    for(i=0; i<nA-nxy; i++) A[i+nxy+i*nA]=e.offDAZn[i];
  }
}

/* --- Private Methods ----------------------------- */
// --- C++ requirement to call Fortran subroutine ----
extern "C"
{
  void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, 
		       int *info);
  void dgetrs_(char *tran, int *n, int *nRHS, double *a, int *lda, 
	       int *ipiv, double *b, int *ldb, int *info);
};

// --- Wrappers to Lapack Fortran subroutines ------
// --- full matrix LU decomposition --------------
int LapackFullSolver::DGETRF_C(int m, int n, double *a, int lda, int *ipiv)
{
  int info = -1;
  dgetrf_(&m, &n, a, &lda, ipiv, &info);
  return info;
}

// --- full matrix back solve using LU results -------------
int LapackFullSolver::DGETRS_C(char tran, int n, int nRHS, double *a, int lda,
			       int *ipiv, double *b, int ldb)
{
  int info = -1;
  dgetrs_(&tran, &n, &nRHS, a, &lda, ipiv, b, &ldb, &info);
  return info;
}

