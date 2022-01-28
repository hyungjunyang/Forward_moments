/*
 * File: Eqn.cc
 * ----------------------------------
 * Implementation for Eqn (Equation) class
 */
#include "Eqn.h"

/* --- Public Methods -------------------------------------- */
// --- Constructors and destructor -------------------
Eqn::Eqn( Grid &grid, Region *regn, Perm *perm, Fluid *fluid ) 
: g(grid), regnPtr(regn), permPtr(perm), fluidPtr(fluid)
{
  // allocate A's elements in subclasses
}

Eqn::Eqn(const Eqn &e) 
: g(e.g), regnPtr(e.regnPtr), permPtr(e.permPtr), fluidPtr(e.fluidPtr)
{
  nx = e.nx; ny = e.ny; 
  nz = e.nz; nNode = e.nNode;
  
  nxPerm = e.nxPerm; nyPerm = e.nyPerm; 
  nzPerm = e.nzPerm; nNodePerm = e.nNodePerm;
  
  unit = e.unit; nA = e.nA; nW = e.nW;

  // --- don't own A, so point to other Eqn(P0Eqn)'s A (all elements) - 
  diagA   = e.diagA; 
  offDAXp = e.offDAXp; offDAYp = e.offDAYp; offDAZp = e.offDAZp;
  offDAXn = e.offDAXn; offDAYn = e.offDAYn; offDAZn = e.offDAZn;
  if(nW > 0) { 
     WW = e.WW;  
     WR = e.WR;  
     RW = e.RW; 
     wEqnInd = e.wEqnInd; 
  } 

  initialize(0.0);  // allocate owned space and assign it to 0
}

Eqn::~Eqn() {
  delete[] RHS; 
  delete[] X; 
  if(fluidPtr->isCompressible()) delete[] X_old; 
}

// --- output A, B, X ------------------------------
ostream &operator<< (ostream &os, Eqn &e) {
  cout << "------X------------\n";
  double aaa = (e.unit == FIELD)? psia_pa : 1;
  for(int k=0; k<e.nz; k++) {
      for(int j=0; j<e.ny; j++) {
          for(int i=0; i<e.nx; i++) {
              cout << i+1 << "  " 
		   << j+1 << "  " 
		   << k+1 << "  " 
	           << e.X[i+(j+k*e.ny)*e.nx]/aaa << endl;
	  }
      }
  }
  return os;
}

/* --- Protected Methods ----------------------------- */
// --- initialize RHS, X, X_old, and assign initial value to X ---
void Eqn::initialize(double X_init)
{
  RHS = new double[nA];
  X   = new double[nA];
  // X = RHS;  // (AAAAA)

  if(fluidPtr->isCompressible()) {    // time dependent
    X_old = new double[nA];
    for(int i=0; i<nNode; i++) X[i] = X_init;
  }
}

// --- calculate the RHS for boundary -----------------------
// --- this is suitable for all high order moments, where the constant in
// --- all B.C is always zero, P0Eqn need to override this method --- 
void Eqn::setBoundaries(bool recalA, const BType *bType, const double *bCond)
{
  int i, j, k;
  int i0, j0, k0, i1, j1, k1; // starting and ending index of each boundary
 
  for(int sur=0; sur<6; sur++) {
    // --- a little tricky, and it is right for all boundarys --- 
    i0 = (sur==1)*(nx-1);  j0 = (sur==3)*(ny-1); k0 = (sur==5)*(nz-1); 
    i1 = (sur!=0)*(nx-1);  j1 = (sur!=2)*(ny-1); k1 = (sur!=4)*(nz-1); 
   
    switch(bType[sur]) {
    case CONST_PRES:   // constant pressure 
      for(k=k0; k<=k1; k++) {
          for(j=j0; j<=j1; j++) {
              for(i=i0; i<=i1; i++) {
	          if(g.getGridType() == POINT_DIS) 
	             RHS[i+(j+k*ny)*nx] = 0.0;
	      }
	  }
      }          // no change to RHS for cell centered grid
      break;
    case CONST_RATE:  // uniform flux
      break;       // no change to RHS for both grids
    default:
      cerr << "UN-Supported boundary condition type: "
	   <<bType[sur]<<endl;
    } 
  }
}

