/*
 * File: FullLUSolver.cc
 * ----------------------------------
 * Implementation for FullLUSolver class
 */
#include "FullLUSolver.h"
#include <math.h>
#include <time.h>
#include "Eqn.h"

/* --- Public Methods ----------------------------------- */
// --- Constructors and destructor ------------
FullLUSolver::FullLUSolver(const Eqn &e) : Solver(e)
{
  A = new double[nA*nA];
  L = new double[nA*nA];
  U = new double[nA*nA];
}

FullLUSolver::~FullLUSolver()
{
  delete[] L; delete[] U; delete[] A;
}

// --- solve AX=B ---------------------------
void FullLUSolver::solve(bool doLU, const Eqn &e, int nRHS, double *RHS)
{
  double startT, endT;
  if(doLU) { 
    // --- A is just assigned or recalculated, need do LU Decomp to A - 
    cout << "NOW, DO LU DECOMP TO MATRIX A --------------- " << endl;
    setupA(e);     // first setup A from Eqn

    // --- do LU Decomposition -----
    startT = clock();
    djlu(A, nA, L, U);
    endT = clock();
    cerr << "LUDecomp time: "<< (endT-startT)/CLOCKS_PER_SEC << " sec" << endl;
 
  }

  // --- backSolves -----------
  backSolveL(L, nRHS, nA, RHS);     // LY=B (Y, B share same space)
  backSolveU(U, nRHS, nA, RHS);     // UX=Y (X, Y share same space)
}

/* --- Protected Methods -------------------------------- */
// --- setup matrix A -----------
void FullLUSolver::setupA(const Eqn &e)
{
  int i;      // Note: A is row wise in C++
  for(i=0; i<e.nA*e.nA; i++) A[i]=0.0;        // all zero
    
  for(i=0; i<nA; i++) A[i+i*nA]=e.diagA[i];   // assign diagonal
  if(nx>1) {                                  // assign off-D X+ and X-
    for(i=1; i<nA; i++) A[i+(i-1)*nA]=e.offDAXp[i-1];
    for(i=0; i<nA-1; i++) A[i+(i+1)*nA]=e.offDAXn[i];
  }
  if(ny>1) {                                  // assign off-D Y+ and Y-
    for(i=nx; i<nA; i++) A[i+(i-nx)*nA]=e.offDAYp[i-e.nx];
    for(i=0; i<nA-nx; i++) A[i+(i+nx)*nA]=e.offDAYn[i];
  }
  if(nz>1) {                                  // assign off-D Z+ and Z-
    int nxy = nx*ny;
    for(i=nxy; i<nA; i++) A[i+(i-nxy)*nA]=e.offDAZp[i-nxy];
    for(i=0; i<nA-nxy; i++) A[i+(i+nxy)*nA]=e.offDAZn[i];
  }
}

/* --- Private Methods ----------------------------- */
// --- do LU decomposition ----------------------------
// --- Note: this method is copied from a book, used as a black box 
// --- so no comment is given, sorry for that ---------------------
int FullLUSolver::djlu(double *a, int n, double *l, double *u)
{
  int i, j, k, w, v, ll;
  
  for(k=0; k<n-1; k++) {
    ll=k*n+k;
    if(fabs(a[ll])+1.0==1.0) { cerr << "Fail\n"; return 0; }
    for(i=k+1; i<n; i++) {
      w=i*n+k; a[w]=a[w]/a[ll];
    }
    for(i=k+1; i<n; i++) {
      w=i*n+k;
      for(j=k+1; j<n; j++) {
	v=i*n+j;
	a[v]=a[v]-a[w]*a[k*n+j];
      }
    }
  }

  for(i=0; i<n; i++) {
    for(j=0; j<i; j++) {
      w=i*n+j; l[w]=a[w]; u[w]=0.0;
    }
    w=i*n+i;
    l[w]=1.0; u[w]=a[w];
    for(j=i+1; j<n; j++) {
      w=i*n+j; l[w]=0.0; u[w]=a[w];
    }
  }
  
  return 1;
}

// --- back solve for L (LX=Y, L is unit lower triangular matrix) --------    
void FullLUSolver::backSolveL(const double *l, int nRHS, 
			    int n, double *x)
{
  int i, j, k;
  double *XX;
  for (k=0; k<nRHS; k++) {
    XX = x+k*n;     // non-blocked, so only deal with one RHS at a time

    //XX[0]/=l[0];
    for(i=1; i<n; i++) {
      for(j=0; j<i; j++) XX[i] -= XX[j]*l[j+i*n];
      //XX[i] /= l[i*n+i];
    }
  }
}

// --- back solve for U (UX=Y, U is upper triangular matrix) --------------
void FullLUSolver::backSolveU(const double *u, int nRHS, int n, double *x)
{
  int i, j, k;
  double *XX;

  for(k=0; k<nRHS; k++) {
    XX = x+k*n;      // non-blocked, so only deal with one RHS at a time
    
    XX[n-1] /= u[n-1+(n-1)*n];
    for(i=n-2; i>=0; i--) {
      for(j=n-1; j>i; j--) XX[i] -= XX[j]*u[j+i*n];
      XX[i] /= u[i*n+i];
    }
  }
}

