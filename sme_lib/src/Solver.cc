/*
 * File: Solver.cc
 * ----------------------------------
 * Implementation for Solver class
 */
#include "Solver.h"
#include <math.h>
#include "Eqn.h"

/* --- Public Methods ---------------------------------------- */
// --- Constructors and destructor ---------------------------
Solver::Solver(const Eqn &e)
{
  debug = false;
  
  nx = e.nx; ny = e.ny; nz = e.nz; nNode = e.nNode;
  nW = e.nW; nA = e.nA;
  nxPerm = e.nxPerm; nyPerm = e.nyPerm; 
  nzPerm = e.nzPerm; nNodePerm = e.nNodePerm; 
  
  if(debug) {
     cout << "In Solver::Solver():" << endl;
     cout << "nxPerm = "    << nxPerm << endl;
     cout << "nyPerm = "    << nyPerm << endl; 
     cout << "nzPerm = "    << nzPerm << endl;
     cout << "nNodePerm = " << nNodePerm  << endl;
  }
}

Solver::~Solver()
{
}

