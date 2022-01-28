 /* File: P0Eqn.cc
 * ----------------------------------
 * Implementation for P0Eqn class
 * Warning:
 * 1) the Strength of wells here is defined as
 *            RHS[ijk] -= w[iw].wStr/g.bdx[ii] * g.bdy[jj] * g.bdz[kk];
 *
 * 2) the matrix A 
 */

#include "P0Eqn.h"

//===== ( Constructor ) =====
P0Eqn::P0Eqn( Grid &grid, Region *regn, Perm *perm, Fluid *fluid,
              int nWells, const Well *w, const Control &c, Unit u, 
              bool debug_) 
: Eqn(grid, regn, perm, fluid), debug(debug_) {
	
  
  if(debug) cout << "P0Eqn::P0Eqn()\n"<<endl;
  
  unit = u; 
  
  nDt    = c.getNDt();
  
  nx        = g.getNx(); 
  ny        = g.getNy(); 
  nz        = g.getNz(); 
  nNode     = g.getnNode();
  nxPerm    = g.getPermNx();
  nyPerm    = g.getPermNy();
  nzPerm    = g.getPermNz(); 
  nNodePerm = g.getNumPermNode();

  nW = 0;
  wEqnInd = new int[nWells];
  for(int i = 0; i < nWells; i++)   // create index mapping (wEqnInd)
      if(w[i].wConType == RATE && w[i].wLens > 1) 
         wEqnInd[nW++] = i;
  nA = nNode + nW;
  if(debug) cout << "nA = " << nA << endl;
  initialize( c.getP_Init() );  // initialize and assign init P to P0
} 
        
//===== ( Destructor ) =====
P0Eqn::~P0Eqn() {
  if(debug) cout << "P0Eqn::~P0Eqn()\n"<<endl;
  delete[] diagA;
  delete[] offDAXp; delete[] offDAXn; 
  delete[] offDAYp; delete[] offDAYn; 
  delete[] offDAZp; delete[] offDAZn; 
  delete[] wEqnInd;
  delete[] WW; delete[] WR; delete[] RW; 
      
  for(int j = 0; j < nDt; j++) {
      if(nx > 1) delete[] dP0dXi[0][j];
      if(ny > 1) delete[] dP0dXi[1][j];
      if(nz > 1) delete[] dP0dXi[2][j];
      if(fluidPtr->isCompressible() ) { 
         delete[] dP0dt[j]; 
         delete[] P0[j]; 
      }
  }   
       
  for(int i = 0; i < 3; i++) 
      delete[] dP0dXi[i]; 
       
  delete[] dP0dt; delete[] P0;
}        
     
//===== ( INITIALIAZE ) =====
void P0Eqn::initialize( double X_init ) {
  Eqn::initialize(X_init);         // initialize inherited members
  initializeA();                   // initialize elements of A
            
  // --- initialize P0's derivatives ---
  // NULL is good for delete ...

  dP0dXi[0] = ( nx > 1 )? new double*[nDt] : NULL;
  if(nx > 1) 
     for(int i = 0; i < nDt; i++) 
         dP0dXi[0][i] = new double[nNode];
                       
  dP0dXi[1] = ( ny > 1 )? new double*[nDt] : NULL;
  if(ny > 1) 
     for(int i = 0; i < nDt; i++) 
         dP0dXi[1][i] = new double[nNode];
  
  dP0dXi[2] = ( nz > 1 )? new double*[nDt] : NULL;
  if(nz > 1) 
     for(int i = 0; i < nDt; i++) 
         dP0dXi[2][i] = new double[nNode];
                            
  bool com = fluidPtr->isCompressible();
  dP0dt = (com)? new double*[nDt] : NULL;
  if(com ) 
     for(int i = 0; i < nDt; i++) 
         dP0dt[i] = new double[nNode];
  
  P0 =(com)? new double*[nDt] : NULL;
  if(com ) 
     for(int i = 0; i < nDt; i++) 
         P0[i] = new double[nNode];
}     

// --- Initialize A ----------------------
void P0Eqn::initializeA() {
  // reservoir diagonal part
  diagA     = new double[nNode];
  
  // reservoir off-diagonal parts
  offDAXp =(nx > 1)? new double[nNode - 1      ] : NULL;
  offDAXn =(nx > 1)? new double[nNode - 1      ] : NULL;
  offDAYp =(ny > 1)? new double[nNode - nx     ] : NULL;
  offDAYn =(ny > 1)? new double[nNode - nx     ] : NULL;
  offDAZp =(nz > 1)? new double[nNode - nx * ny] : NULL;
  offDAZn =(nz > 1)? new double[nNode - nx * ny] : NULL;
  // well parts
  WW =(nW > 0)? new double[nW*nW]    : NULL;
  RW =(nW > 0)? new double[nW*nNode] : NULL;
  WR =(nW > 0)? new double[nNode*nW] : NULL;
}

//===== ( Solve ) =====
void P0Eqn::solve(int nWells, const Well *w, const Control &c, 
                  Solver &s, int iter) {

  for(int i = 0; i < nNode; ++i) RHS[i] = 0.;

  // --- assign X_old (if time dependent) ---
  // if incompressible, then no time dependent
  if( fluidPtr->isCompressible() ) {
      for(int i = 0; i < nNode; i++ ) 
          X_old[i] = X[i];
  }

  // --- setup A and B ---
  bool recalA = ( (iter == 0) || 
                  (c.getDt( iter ) != c.getDt( iter-1 )) || 
                  fluidPtr->isCompressible() 
		);
             
  // A for internal points
  if( recalA  ) {
      if(debug) {
         cout << endl;
         cout << "calcA() in P0Eqn.cc" << endl << endl;
         cout << "recalA = " <<  recalA << endl;
      }
      calcA();
  }
  
  // B for internal points
  // calcRHS( NULL ); 
       
  // boundary points
  setBoundaries( recalA, c.getBType(), c.getBCond() );
               
  // well blocks
  setWells( recalA, nWells, w );
              
  // time-derivatives 
  if(fluidPtr->isCompressible() ) 
     calcAccu( recalA, c.getDt(iter) );
               
  // --- Solve for P0 ---
  // pass to solver, solve it.
  
 /*
  ofstream os("diagA1", ios::out);
  for(int i = 0; i < nNode; i++) {
      os << diagA[i] << ' '<< RHS[i] << endl;
  }
  os << endl;
  exit(0);
  
  for(int i = 0; i < nNode - 1; i++) 
      os << offDAXp[i] << ' '<< offDAXn[i] << endl;
  os << endl;
	  
  for(int i = 0; i < nNode - nx; i++) 
      os << offDAYp[i] <<' ' << offDAYn[i] << endl;
  os << endl;
	  
  for(int i = 0; i < nNode - nx*ny; i++) 
      os << offDAZp[i] <<' ' << offDAZn[i] << endl;
  os << endl;

  exit(0);  
  */
  
  s.solve( recalA, *this, 1, RHS );
  
  /*
  ofstream os("p0_1p.out", ios::out);
  for(int i = 0; i < nNode; ++i) {
      os<< RHS[i]<< endl;
  }
  os<<endl;
  os.close();
  */
  for(int i = 0; i < nNode; i++) {
      X[i] = RHS[i];
  }

  // --- post-processing and updating -------------
  // record P0 for each timestep
  if(fluidPtr->isCompressible() ) {
     for(int i = 0; i < nNode; i++) 
         P0[iter][i] = X[i];
     // update den, poro...
     fluidPtr->update( X, false );
  }             
         
  // --- calculate P0's derivatives -------------
  // P0's derivatives
  calcDerivatives(iter, c.getDt(iter));
  // calcDerivatives(iter, c.dt[iter], c.bType, c.bCond);
}            

//===== ( Solve ) =====
void P0Eqn::solve(int nWells, const Well *w, const Control &c, 
                  Solver &s, int iter, int num_meas, 
		  int* i_pmeas, int* j_pmeas, int* k_pmeas, double *p_pmeas) {

  for(int i = 0; i < nNode; ++i) RHS[i] = 0.;

  // --- assign X_old (if time dependent) ---
  if( fluidPtr->isCompressible() ) {
      for(int i = 0; i < nNode; i++ ) 
          X_old[i] = X[i];
  }

  // --- setup A and B ---
  bool recalA = ( (iter == 0) || 
                  (c.getDt( iter ) != c.getDt( iter-1 )) || 
                  fluidPtr->isCompressible() 
		);
             
  // A for internal points
  calcA();
  
  // boundary points
  setBoundaries( true, c.getBType(), c.getBCond() );
               
  // well blocks
  setWells( true, nWells, w );
 
  // Conditioning part 
  for(int m = 0; m < num_meas; m++ ) {
      int i   = i_pmeas[m];
      int j   = j_pmeas[m];
      int k   = k_pmeas[m];
      int ijk = g.getIndex(i, j, k);
      double vol = g.getBdx( i ) * g.getBdy( j ) * g.getBdz( k );
      //cout << "vol"  << vol << endl;
      //diagA[ijk] -= 1.0e14/vol;
      //  RHS[ijk] -= 1.0e14  *  p_pmeas[m] / vol;
      diagA[ijk] = -1.0e14/vol;
        RHS[ijk] = -1.0e14  *  p_pmeas[m] / vol;
        //cout <<"cond:" << i << ' '<<j <<' '<< k<<' '<<ijk<<' '
        //     << diagA[ijk]<<' '<<RHS[ijk]<<endl;
  }
  
  // time-derivatives 
  if(fluidPtr->isCompressible() ) 
     calcAccu( recalA, c.getDt(iter) );
  /* 
  ofstream os("diagA", ios::out);
  for(int i = 0; i < nNode; i++) {
      os << diagA[i] << ' '<< RHS[i] << endl;
  }
  for(int i = 0; i < nNode - 1; i++) 
      os << offDAXp[i] << ' '<< offDAXn[i] << endl;
  for(int i = 0; i < nNode - nx; i++) 
      os << offDAYp[i] <<' ' << offDAYn[i] << endl;
  exit(0); 
  */ 
  // --- Solve for P0 ---
  /*
  cout << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(6);
  cout << endl;
  for(int i = 0; i < nNode; i++) {
      cout << offDAXn[i - 1 ]  << ' '
           << offDAYn[i - nx]  << ' '
           <<   diagA[i]  << ' '
           << offDAXp[i]  << ' '
           << offDAYp[i]  ;
      //     <<     RHS[i];
      cout << endl;
  }
  */
  s.solve( true, *this, 1, RHS );

  for(int i = 0; i < nNode; i++) 
      X[i] = RHS[i];
                     
  // --- post-processing and updating -------------
  if(fluidPtr->isCompressible() ) {
     for(int i = 0; i < nNode; i++) 
         P0[iter][i] = X[i];
     fluidPtr->update( X, false );
  }             
         
  calcDerivatives(iter, c.getDt(iter));
}            

//=========Recalculate A ==========================
void P0Eqn::recalcA(const BType *bType, int nWells, const Well *w, double dt) {
  cout << "RE-CALCULATE A IN POEQN ----------------------\n";
  cout << "ent at recalcA() of P0EQn\n";
  exit(0);
  int i, ijk;  double vol;
         
  // --- for internal points ---
  calcA();
           
  // --- for boundary points ---
  /* in order to correctly deal with shared boundary points (corners),
   * we need to separate two type of boundaries, const pressure boundary 
   * has the priority
   */        
  // --- first deal with uniform flux boundaries     ---
  if( nx > 1 && bType[0] == CONST_RATE ) aBoundXp( bType[0] ); 
  if( nx > 1 && bType[1] == CONST_RATE ) aBoundXn( bType[1] ); 
  if( ny > 1 && bType[2] == CONST_RATE ) aBoundYp( bType[2] ); 
  if( ny > 1 && bType[3] == CONST_RATE ) aBoundYn( bType[3] ); 
  if( nz > 1 && bType[4] == CONST_RATE ) aBoundZp( bType[4] ); 
  if( nz > 1 && bType[5] == CONST_RATE ) aBoundZn( bType[5] ); 
                                                   
  // --- then deal with constant pressure boundaries  ---    
  if( nx > 1 && bType[0] == CONST_PRES ) aBoundXp( bType[0] ); 
  if( nx > 1 && bType[1] == CONST_PRES ) aBoundXn( bType[1] ); 
  if( ny > 1 && bType[2] == CONST_PRES ) aBoundYp( bType[2] ); 
  if( ny > 1 && bType[3] == CONST_PRES ) aBoundYn( bType[3] ); 
  if( nz > 1 && bType[4] == CONST_PRES ) aBoundZp( bType[4] ); 
  if( nz > 1 && bType[5] == CONST_PRES ) aBoundZn( bType[5] ); 

  // --- for wells ---------------
  for( int iw = 0; iw < nWells; iw++) {
       switch( w[iw].wConType ) {
          case RATE:        // constant Rate  (no change)
               break;
          case PRES:        // contant pressure
              for(i = 0; i < w[iw].wLens; i++) {
                  ijk =    w[iw].wLocIJK[0][i] 
                     + (   w[iw].wLocIJK[1][i] 
                         + w[iw].wLocIJK[2][i]*ny
                       )*nx;
                  vol =  g.getBdx( w[iw].wLocIJK[0][i] )
                       * g.getBdy( w[iw].wLocIJK[1][i] )
                       * g.getBdz( w[iw].wLocIJK[2][i] );  
                  diagA[ijk] -=  w[iw].wTrans[i]/vol;   
              }
          break;
       }
  }
  // --- for time-dependent part ---
  if( fluidPtr->isCompressible() )
      for(i = 0; i < nNode; i++) 
          diagA[i] -= fluidPtr->getPoro_crf( i )/dt;
}         
      
//============ Calculate Matrix A ===============
void P0Eqn::calcA() {
  int i, j, k, ijk;
  double tmp_p, tmp_n, dd;
  double perm_interface[3];
                
  // --- set A to zero ----------
  for(i = 0; i < nNode; i++) 
      diagA[i]   = 0.0;
  if( nx > 1 ) 
      for(i = 0; i < nNode - 1; i++) 
          offDAXp[i] = offDAXn[i] = 0.0;
  if( ny > 1 ) 
      for(i = 0; i < nNode - nx; i++) 
          offDAYp[i] = offDAYn[i] = 0.0;
  if( nz > 1 ) 
      for(i = 0; i < nNode - nx*ny; i++) 
          offDAZp[i] = offDAZn[i] = 0.0;
                   
  // --- X direction ----------------------
  if( nx > 1 ) {
      for(k = 0; k < nz; k++ ) {
          for(j = 0; j < ny; j++ ) {
              for(i = 1; i < nx - 1; i++ ) {
                  ijk = g.getIndex(i, j, k );
                  permPtr->Perm_x( i, j, k, perm_interface );
                  offDAXp[ijk  ] = perm_interface[2] / g.getBdx( i )
                                                     / g.getDx(i   );
                  offDAXn[ijk-1] = perm_interface[0] / g.getBdx( i )
                                                     / g.getDx(i - 1);
                    diagA[ijk] -= ( offDAXp[ijk] + offDAXn[ijk-1] );
              }      
          }            
      }                  
  }                    
  // --- Y direction ----------------------
  if( ny > 1) {              
      for(j = 1; j < ny - 1; j++) {
          for(k = 0; k < nz; k++) {
              for(i = 0; i < nx; i++) {
                  ijk = g.getIndex( i, j, k ); 
                  permPtr->Perm_y( i, j, k, perm_interface );
                  offDAYp[ijk   ] = perm_interface[2] / g.getBdy( j )
                                                      / g.getDy( j  );
                  offDAYn[ijk-nx] = perm_interface[0] / g.getBdy( j )
                                                      / g.getDy(j - 1);
                  diagA[ijk] -= ( offDAYp[ijk] + offDAYn[ijk-nx] );
              }                               
          }                                       
      }
  }                                        
                                                  
  // --- Z direction (include gravity)--------------
  if( nz > 1 ) {                      
      for(k = 1; k < nz - 1; k++) {
          for(j = 0; j < ny; j++) {
              for(i = 0; i < nx; i++) {
                  ijk = g.getIndex( i, j, k ); 
                  permPtr->Perm_z( i, j, k, perm_interface );
                  offDAZp[ijk      ] = perm_interface[2]/g.getBdz(k)
                                                        /g.getDz(k);
                  offDAZn[ijk-nx*ny] = perm_interface[0]/g.getBdz(k)
                                                        /g.getDz(k-1);
                  diagA[ijk] -= ( offDAZp[ijk] + offDAZn[ijk-nx*ny] );

                  //offDAZp[ijk      ] += perm_interface[1]*fluidPtr->getDen_cf(ijk)
                  //                      / 2.0 /g.getBdz(k);
                  //offDAZn[ijk-nx*ny] += -perm_interface[1]*fluidPtr->getDen_cf(ijk)
                  //                      / 2.0 /g.getBdz(k);
                  RHS[ijk] = -fluidPtr->getGamma(ijk) *
                            ( perm_interface[2] - perm_interface[0] ) 
                            / g.getBdz(k);
              }                     
          }                       
      }                                      
  }                                         
} 

// aBound means the boundary of Matrix A
//===== ( i == 0 ) =====
void P0Eqn::aBoundXp( BType bType ) {
  int j, k, ijk;
  double gdx  = g.getDx( 0),  gdx2 = gdx  * gdx;
  double gbdx = g.getBdx( 0), gbdx2 = gbdx * gbdx;
  double tmp;
  double perm_interface[3];
  if( g.getGridType() == POINT_DIS ) {
      switch( bType ) {
        case CONST_PRES: // contant pressure
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(0, j, k); // (j+k*ny)*nx;
                    if( j < ny - 1 ) offDAYp[ijk      ] = 0.0;
                    if( j > 0)       offDAYn[ijk - nx ] = 0.0;
                    if( k < nz - 1 ) offDAZp[ijk      ] = 0.0;
                    if( k > 0 )      offDAZn[ijk-nx*ny] = 0.0;
                                       diagA[ijk      ] = 1.0;
               }
          }
          break;
        case CONST_RATE: // uniform flux (No Flow Boundary)
          tmp = 2./gdx2;     
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(0, j, k); 
                    permPtr->Perm_x( 0, j, k, perm_interface );
                    offDAXp[ijk]  = tmp * perm_interface[2];
                      diagA[ijk] -= offDAXp[ijk];
               }
          }
          break;
      } 
  }       
  else { // cell centered
      tmp = 1.0/gbdx/gdx;
      switch( bType ) {
        case CONST_PRES:                       // constant pressure
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(0, j, k);
                    permPtr->Perm_x( 0, j, k, perm_interface );
                    offDAXp[ijk]  = tmp * perm_interface[2]; 
                      diagA[ijk] -= offDAXp[ijk] + 2. * perm_interface[0]/gbdx2; 
               }
          }
          break;
        case CONST_RATE:                       // constant flux rate     
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) { 
                    ijk = g.getIndex(0, j, k);       
                    permPtr->Perm_x( 0, j, k, perm_interface );
                    offDAXp[ijk]   = tmp * perm_interface[2];
                      diagA[ijk]  -= offDAXp[ijk];
               }
          }
          break;
      }
  }
}

//===== ( i == nx -1 ) =====
void P0Eqn::aBoundXn( BType bType ) {
  int j, k, ijk;
  double gdx  = g.getDx( nx - 2 ),  gdx2 = gdx  * gdx;
  double gbdx = g.getBdx( nx - 1 ), gbdx2 = gbdx * gbdx;
  double tmp;
  double perm_interface[3];
  if( g.getGridType() == POINT_DIS ) {
      switch( bType ) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(nx - 1, j, k); 
                    if( j < ny - 1) offDAYp[ijk]       = 0.0;
                    if( j > 0)      offDAYn[ijk-nx]    = 0.0;
                    if( k < nz - 1) offDAZp[ijk]       = 0.0;
                    if( k > 0)      offDAZn[ijk-nx*ny] = 0.0;
                                      diagA[ ijk ]     = 1.0;
               }
          }
          break;
        case CONST_RATE:
          tmp = 2./gdx2;
          for(k = 0; k < nz; k++) {
              for(j = 0; j < ny; j++) {
                  ijk = g.getIndex(nx - 1, j, k); 
                  permPtr->Perm_x( nx - 1, j, k, perm_interface );
                  offDAXn[ijk-1] = tmp * perm_interface[0];
                    diagA[ijk]  -= offDAXn[ijk-1];
              }
          }
          break;
      }
  }
  else {
      tmp = 1.0/gbdx/gdx;
      switch(bType) {
        case CONST_PRES:                  
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) { 
                    ijk = g.getIndex(nx - 1, j, k); 
                    permPtr->Perm_x( nx - 1, j, k, perm_interface );
                    offDAXn[ijk-1] = tmp * perm_interface[0];
                       diagA[ijk] -= offDAXn[ijk-1] + 2. * perm_interface[2]/gbdx2;
               }
          }
          break;
        case CONST_RATE:                  
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) { 
                    ijk = g.getIndex(nx - 1, j, k); 
                    permPtr->Perm_x( nx - 1, j, k, perm_interface );  
                    offDAXn[ijk-1] = tmp * perm_interface[0];
                      diagA[ijk ] -= offDAXn[ijk-1]; 
               }
          }
        break;
      }
  } 
}

//===== ( j == 0 ) =====
void P0Eqn::aBoundYp(BType bType) {
  int i, k, ijk;
  double gdy  = g.getDy( 0 ),  gdy2 = gdy  * gdy;
  double gbdy = g.getBdy( 0 ), gbdy2 = gbdy * gbdy;
  double tmp;
  double perm_interface[3];
  if( g.getGridType() == POINT_DIS ) {
      switch( bType ) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, 0, k);
                    if(i < nx-1) offDAXp[ijk      ] = 0.0;
                    if(i > 0   ) offDAXn[ijk-1    ] = 0.0;
                    if(k < nz-1) offDAZp[ijk      ] = 0.0;
                    if(k > 0   ) offDAZn[ijk-nx*ny] = 0.0;
                                   diagA[ijk      ] = 1.0;
               }
          }
          break;
        case CONST_RATE:
          tmp = 2./gdy2;
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, 0, k);
                    permPtr->Perm_y( i, 0, k, perm_interface );
                    offDAYp[ijk]  = tmp * perm_interface[2];
                      diagA[ijk] -= offDAYp[ijk]; 
               }
          }
          break;
      }
  }
  else {
      tmp = 1.0/gbdy/gdy;
      switch( bType ) {
        case CONST_PRES:                  
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) { 
                    ijk = g.getIndex(i, 0, k);
                    permPtr->Perm_y( i, 0, k, perm_interface ); 
                    offDAYp[ijk]  = tmp * perm_interface[2];
                      diagA[ijk] -= offDAYp[ijk] + 2. * perm_interface[0]/gbdy2;
               }
          }
          break;
        case CONST_RATE:                  
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) { 
                    ijk = g.getIndex(i, 0, k);
                    permPtr->Perm_y( i, 0, k, perm_interface );  
                    offDAYp[ijk]  = tmp * perm_interface[2];
                      diagA[ijk] -= offDAYp[ijk]; 
               }
          }
          break;
      }
  }
}

//===== ( j == ny - 1) =====
void P0Eqn::aBoundYn( BType bType ) {
  int i, k, ijk;
  double gdy  = g.getDy( ny - 2 ),  gdy2 = gdy  * gdy;
  double gbdy = g.getBdy( ny - 1 ), gbdy2 = gbdy * gbdy;
  double tmp;
  double perm_interface[3];  
  if( g.getGridType() == POINT_DIS) {
      switch( bType ) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, ny - 1, k);
                    if( i < nx - 1) offDAXp[ijk]   = 0.0;
                    if( i > 0     ) offDAXn[ijk-1] = 0.0;
                    if( k < nz - 1) offDAZp[ijk]   = 0.0;
                    if( k > 0     ) offDAZn[ijk-1] = 0.0;
                                      diagA[ijk]   = 1.0;
               }
          }
          break;
        case CONST_RATE:
          tmp = 2./gdy2;
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, ny - 1, k);
                    permPtr->Perm_y( i, ny - 1, k, perm_interface );      
                    offDAYn[ijk-nx] = tmp * perm_interface[0];
                      diagA[ijk]   -= offDAYn[ijk-nx]; 
               }
          }
          break;
      }      
  }       
  else { 
      tmp = 1.0/gbdy/gdy;
      switch(bType) {
        case CONST_PRES:        
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) { 
                    ijk = g.getIndex(i, ny - 1, k); 
                    permPtr->Perm_y( i, ny - 1, k, perm_interface );
                    offDAYn[ijk-nx] = tmp * perm_interface[0];
                      diagA[ijk  ] -= offDAYn[ijk-nx] + 2. * perm_interface[2]/gbdy2;
               }
          }
          break;
        case CONST_RATE:      
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) { 
                    ijk = g.getIndex(i, ny - 1, k); 
                    permPtr->Perm_y( i, ny - 1, k, perm_interface );
                    offDAYn[ijk-nx] = tmp * perm_interface[0];
                      diagA[ijk  ] -= offDAYn[ijk-nx]; 
               }
          } 
          break;
      }
  }
}

//===== ( k == 0) =====
void P0Eqn::aBoundZp( BType bType ) {
  int i, j, ijk;
  double gdz  = g.getDz( 0 ),  gdz2 = gdz  * gdz;
  double gbdz = g.getBdz( 0 ), gbdz2 = gbdz * gbdz;
  double tmp;
  double perm_interface[3];
  if( g.getGridType() == POINT_DIS ) {
      switch( bType ) {
        case CONST_PRES:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, 0); 
                    if(i < nx-1) offDAXp[ijk   ] = 0.0;
                    if(i > 0   ) offDAXn[ijk-1 ] = 0.0;
                    if(j < ny-1) offDAYp[ijk   ] = 0.0;
                    if(j > 0   ) offDAYn[ijk-nx] = 0.0;
                                   diagA[ijk   ] = 1.0; 
               }
          }
          break;
        case CONST_RATE:
          tmp = 2./gdz2;     
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, 0);
                    permPtr->Perm_z( i, j, 0, perm_interface );
                    offDAZp[ijk]  = tmp * perm_interface[2];
                      diagA[ijk] -= offDAZp[ijk];
               }
          }
          break;
      }    
  }          
  else {       
      tmp = 1.0/gbdz/gdz;
      switch(bType) {
        case CONST_PRES:      
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) { 
                    ijk = g.getIndex(i, j, 0); 
                    permPtr->Perm_z( i, j, 0, perm_interface );
                    offDAZp[ijk]  = tmp * perm_interface[2];
                      diagA[ijk] -= offDAZp[ijk] + 2. * perm_interface[0]/gbdz2;
               }
          }
          break;
        case CONST_RATE:       
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) { 
                    ijk = g.getIndex(i, j, 0); 
                    permPtr->Perm_z( i, j, 0, perm_interface );
                    offDAZp[ijk]  = tmp * perm_interface[2];
                      diagA[ijk] -= offDAZp[ijk];
               }
          }
          break;
      }
  }
}

//===== ( k == nz - 1) =====
void P0Eqn::aBoundZn(BType bType) {
  int i, j, ijk;
  double gdz  = g.getDz( nz - 2 ),  gdz2 = gdz  * gdz;
  double gbdz = g.getBdz( nz - 1 ), gbdz2 = gbdz * gbdz;
  double tmp;
  double perm_interface[3];
  if( g.getGridType() == POINT_DIS ) {
      switch( bType ) {
        case CONST_PRES:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, nz - 1); 
                    if( i < nx-1) offDAXp[ijk   ] = 0.0;
                    if( i > 0   ) offDAXn[ijk-1 ] = 0.0;
                    if( j < ny-1) offDAYp[ijk   ] = 0.0;
                    if( j > 0   ) offDAYn[ijk-nx] = 0.0;
                                    diagA[ijk   ] = 1.0;
               }
          }
          break;
        case CONST_RATE:
          tmp = 2./gdz2;
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, nz - 1);
                    permPtr->Perm_z( i, j, nz - 1, perm_interface );
                    offDAZn[ijk-nx*ny] = tmp * perm_interface[0];
                      diagA[ijk]      -= offDAZn[ijk-nx*ny];
               }
          }
          break;
      }
  }       
  else {
      tmp = 1.0/gbdz/gdz;
      switch(bType) {
        case CONST_PRES:        
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, nz - 1);
                    permPtr->Perm_z( i, j, nz - 1, perm_interface );
                    offDAZn[ijk-nx*ny] = tmp * perm_interface[0];
                      diagA[ijk]      -= offDAZn[ijk-nx*ny] 
                                       + 2. * perm_interface[2]/gbdz2;
               }
          }
          break;
        case CONST_RATE:        
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) { 
                    ijk = g.getIndex(i, j, nz - 1);
                    permPtr->Perm_z( i, j, nz - 1, perm_interface );
                    offDAZn[ijk-nx*ny] = tmp * perm_interface[0];
                      diagA[ijk]      -= offDAZn[ijk-nx*ny]; 
               }
          }
          break;
      }
  } 
}

//===== ( calcRHS(const BType* bType ) =====
void P0Eqn::calcRHS(const BType* bType) {
    cout<<"no implementation!"<<endl;
}

//===== ( Set Boundaries ) =====
void P0Eqn::setBoundaries( bool recalA, const BType *bType,
	       const double *bCond ) {
                           
  // --- first deal with uniform flux boundaries  ---
  if( nx > 1 && bType[0] == CONST_RATE ) {
      if(recalA) aBoundXp( bType[0] ); 
                 bBoundXp( bType[0], bCond[0] ); 
  }
          
  if( nx > 1 && bType[1] == CONST_RATE ) {
      if(recalA) aBoundXn( bType[1] ); 
      bBoundXn( bType[1], bCond[1] ); 
  }
      
  if( ny > 1 && bType[2] == CONST_RATE ) {
      if(recalA) aBoundYp( bType[2] ); 
      bBoundYp( bType[2], bCond[2] ); 
  }
    
  if( ny > 1 && bType[3] == CONST_RATE ) {
      if(recalA) aBoundYn( bType[3] ); 
      bBoundYn( bType[3], bCond[3] ); 
  }
         
  if( nz > 1 && bType[4] == CONST_RATE ) {
      if(recalA) aBoundZp( bType[4] ); 
      bBoundZp( bType[4], bCond[4] ); 
  }
            
  if( nz > 1 && bType[5] == CONST_RATE ) {
      if(recalA) aBoundZn( bType[5] ); 
      bBoundZn( bType[5], bCond[5] ); 
  }
       
  // --- then deal with constant pressure boundaries  ---
  if( nx > 1 && bType[0] == CONST_PRES ) {
      if(recalA) aBoundXp( bType[0] ); 
      bBoundXp( bType[0], bCond[0] ); 
  }
      
  if( nx > 1 && bType[1] == CONST_PRES ) {
      if(recalA) aBoundXn( bType[1] ); 
      bBoundXn( bType[1], bCond[1] ); 
  }

  if( ny > 1 && bType[2] == CONST_PRES ) {
      if(recalA) aBoundYp( bType[2] ); 
      bBoundYp( bType[2], bCond[2] ); 
  }
         
  if( ny > 1 && bType[3] == CONST_PRES ) {
      if(recalA) aBoundYn( bType[3] ); 
      bBoundYn( bType[3], bCond[3] ); 
  }

  if( nz > 1 && bType[4] == CONST_PRES ) {
      if(recalA) aBoundZp( bType[4] ); 
      bBoundZp( bType[4], bCond[4] ); 
  }
       
  if( nz > 1 && bType[5] == CONST_PRES ) {
      if(recalA) aBoundZn( bType[5] ); 
      bBoundZn( bType[5], bCond[5] ); 
  }
}

//===== ( i == 0) =====
void P0Eqn::bBoundXp( BType bType, double bCond ) {
  int j, k, ijk;
  double gdx  = g.getDx( 0 ),  gdx2 = gdx  * gdx;
  double gbdx = g.getBdx( 0 ), gbdx2 = gbdx * gbdx;
  double tmp;
  double perm_interface[3];
  if( g.getGridType() == POINT_DIS ) {
      switch( bType ) {
        case CONST_PRES:  // contant pressure
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex( 0, j, k ); 
                    RHS[ ijk ] = bCond;          
               }
          }
          break; 
        case CONST_RATE: // contant flux 
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(0, j, k ); 
                    permPtr->Perm_x( 0, j, k, perm_interface );
                    RHS[ijk] -= 2.0 * perm_interface[2] * bCond /
                                      ( perm_interface[1] * gdx );
               }
          }
          break;
      }
  }
  else {
      switch(bType) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(0, j, k ); 
                    permPtr->Perm_x( 0, j, k, perm_interface );  
                    RHS[ijk] -= 2.0 * perm_interface[0] * bCond /
                                      gbdx2; 
               }
          }
          break;
        case CONST_RATE:
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(0, j, k ); 
                    RHS[ijk] -= bCond / gbdx;
               }
          }
          break;
      }
  }
}

//===== ( i = nx - 1) =====
void P0Eqn::bBoundXn( BType bType, double bCond ) {
  int j, k, ijk;
  double gdx  = g.getDx( nx - 2 ),  gdx2 = gdx  * gdx;
  double gbdx = g.getBdx( nx - 1 ), gbdx2 = gbdx * gbdx;
  double tmp;
  double perm_interface[3];
  if( g.getGridType() == POINT_DIS ) {
      switch( bType ) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex( nx - 1, j, k );
                    RHS[ ijk ] = bCond;
               }
          }
          break;
        case CONST_RATE:
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(nx - 1, j, k );
                    permPtr->Perm_x( nx - 1, j, k, perm_interface );
                    RHS[ijk] -= 2.0 * perm_interface[0] * bCond /
                                     ( perm_interface[1] * gdx );
               }
          }
          break;
      }   
  }           
  else {     
      switch( bType ) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(nx - 1, j, k ); 
                    permPtr->Perm_x( nx - 1, j, k, perm_interface );
                    RHS[ijk] -= 2.0 * perm_interface[2] * bCond /
                                      gbdx2;
               }
          }
          break;
        case CONST_RATE:
          for( k = 0; k < nz; k++) {
               for( j = 0; j < ny; j++) {
                    ijk = g.getIndex(nx - 1, j, k );
                    RHS[ijk] -= bCond / gbdx;
               } 
          }
          break;
      }
  }
}

//===== ( J = 0 ) =====
void P0Eqn::bBoundYp( BType bType, double bCond ) {
  int i, k, ijk;
  double gdy  = g.getDy( 0 ),  gdy2 = gdy  * gdy;
  double gbdy = g.getBdy( 0 ), gbdy2 = gbdy * gbdy;
  double tmp;
  double perm_interface[3];  
  if( g.getGridType() == POINT_DIS) {
      switch( bType ) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex( i, 0, k );
                    RHS[ ijk ] = bCond;
               }
          }
          break;
        case CONST_RATE:
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, 0, k );
                    permPtr->Perm_y( i, 0, k, perm_interface );
                    RHS[ijk] -= 2.0 * perm_interface[2] * bCond /
                                ( perm_interface[1] * gdy );
               }
          }
          break;
      }
  }     
  else {
      switch( bType ) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
              for( i = 0; i < nx; i++) {
                   ijk = g.getIndex( i, 0, k );
                   permPtr->Perm_y( i, 0, k, perm_interface );
                   RHS[ijk] -= 2.0 * perm_interface[0] * bCond /
                                      gbdy2;
              }
          }
          break;
        case CONST_RATE:
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex( i, 0, k );
                    RHS[ijk] -= bCond / gbdy;
               }
          }
          break;
      }
  }
}

//===== ( j == ny -1) =====
void P0Eqn::bBoundYn( BType bType, double bCond ) {
  int i, k, ijk;
  double gdy  = g.getDy( ny - 2 ),  gdy2 = gdy  * gdy;
  double gbdy = g.getBdy( ny - 1 ), gbdy2 = gbdy * gbdy;
  double tmp;
  double perm_interface[3];  
  if( g.getGridType() == POINT_DIS) {
      switch( bType ) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex( i, ny - 1, k ); 
                    RHS[ ijk ] = bCond;
               }
          }
          break;
        case CONST_RATE:
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, ny - 1, k );  
                    permPtr->Perm_y( i, ny - 1, k, perm_interface );
                    RHS[ijk] -= 2.0 * perm_interface[0] * bCond /
                                ( perm_interface[1] * gdy );
               }
          }
          break;
      }
  }
  else {
      switch( bType ) {
        case CONST_PRES:
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, ny - 1, k );
                    permPtr->Perm_y( i, ny - 1, k, perm_interface );
                    RHS[ijk] -= 2.0 * perm_interface[2] * bCond /
                                      gbdy2;
               }
          }
          break;
        case CONST_RATE:
          for( k = 0; k < nz; k++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex( i, ny - 1, k );
                    RHS[ijk] -= bCond / gbdy;
               }
          }
          break;
      }
  }
}

//===== ( k == 0 ) =====
//a = 9.8 * den, b = 9.8*den*cf
void P0Eqn::bBoundZp( BType bType, double bCond ) {
  int i, j, ijk;
  double ddzp = g.getDz( 0 ), ddzn, ddz;
  double gdz  = g.getDz( 0 ),  gdz2 = gdz  * gdz;
  double gbdz = g.getBdz( 0 ), gbdz2 = gbdz * gbdz;
  double  perm_interface[3];
  if( g.getGridType() == POINT_DIS) {
      switch( bType ) {
        case CONST_PRES:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex( i, j, 0 );
                    RHS[ijk] = bCond;
               }
          }
          break;
        case CONST_RATE:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, 0 );
                    permPtr->Perm_z( i, j, 0, perm_interface );
                    RHS[ijk] -= 2.0 * perm_interface[2] / gdz *
                                   ( bCond/perm_interface[1] + fluidPtr->getGamma(ijk) );
               }
          }
          break;
      }
  }
  else {
      switch( bType ) {
        case CONST_PRES:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, 0 ); 
                    permPtr->Perm_z( i, j, 0, perm_interface );    
                    RHS[ijk] -= ( 2.0 * perm_interface[0] * bCond / gbdz2
                                 + fluidPtr->getGamma(ijk)
                                   * ( perm_interface[2] - perm_interface[0]) /gbdz
                                );    
                        
               }
          }
          break;
        case CONST_RATE:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, 0 );
                    permPtr->Perm_z( i, j, 0, perm_interface );
                    RHS[ijk] -= (  bCond / gbdz
                                 + fluidPtr->getGamma(ijk) * perm_interface[2] /gbdz
                                );
               }
          }
          break;
      }
  }
}

//===== (k == nz - 1 ) =====
void P0Eqn::bBoundZn( BType bType, double bCond ) {
  int i, j, ijk;
  double gdz  = g.getDz( nz - 2 ),  gdz2 = gdz  * gdz;
  double gbdz = g.getBdz( nz - 1 ), gbdz2 = gbdz * gbdz;  
  double perm_interface[3];  
  if( g.getGridType() == POINT_DIS ) {
      switch(bType) {
        case CONST_PRES:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex( i, j, nz - 1 ); 
                    RHS[ijk] = bCond;
               }
          }
          break;
        case CONST_RATE:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, nz - 1 );
                    permPtr->Perm_z( i, j, nz - 1, perm_interface );
                    RHS[ijk] -= 2.0 *  perm_interface[0] *
                       ( bCond/perm_interface[1] - fluidPtr->getGamma(ijk) ) / gdz;
               }
          }
          break;
      }            
  }                  
  else {
      switch(bType) {
        case CONST_PRES:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, nz - 1 ); 
                    permPtr->Perm_z( i, j, nz - 1, perm_interface );    
                    RHS[ijk] -= ( 2.0 * perm_interface[2] * bCond / gbdz2
                                 + fluidPtr->getGamma(ijk)
                                   * ( perm_interface[2] - perm_interface[0]) /gbdz
                                ); 
               }
          }
          break;
        case CONST_RATE:
          for( j = 0; j < ny; j++) {
               for( i = 0; i < nx; i++) {
                    ijk = g.getIndex(i, j, nz - 1 );
                    permPtr->Perm_z( i, j, nz - 1, perm_interface ); 
                    RHS[ijk] -= (  bCond / gbdz
                                 - fluidPtr->getGamma(ijk) * perm_interface[0] /gbdz    
                                );
               }
          }
          break;
      }
  } 
}

//===== ( setWells ) =====
void P0Eqn::setWells( bool recalA, int nWells, const Well *w ) {
  int i, iw, ijk, ii, jj, kk;
  double vol, gra;
  for( iw = 0; iw < nWells; iw++ ) {
       switch( w[iw].wConType ) {
        case RATE:
          // --- only deal with single block penetration now (AAAA) ---
          // for each well block
          for(i = 0; i < w[iw].wLens; i++) {
              //ijk,vol: well block index and volume
              ii  = w[iw].wLocIJK[0][i];
              jj  = w[iw].wLocIJK[1][i];
              kk  = w[iw].wLocIJK[2][i];
              ijk = g.getIndex(ii, jj, kk);
              vol = g.getBdx( ii ) * g.getBdy( jj ) * g.getBdz( kk );
              RHS[ijk] -= w[iw].wStr/vol;
	      if(debug) {
                 cout << "w[iw].wStr = " << w[iw].wStr<<endl;
                 cout << "bdx = " 
                      << g.getBdx( ii ) <<' '
                      << g.getBdy( jj ) <<' ' 
                      << g.getBdz( kk ) <<endl;
                 cout << "the well strength = " << w[iw].wStr/vol << endl;
                 cout << "RHS = " << RHS[ijk] <<endl;
	      }
	      //cout << "the well strength = " << w[iw].wStr<<' '<<vol<< endl;
	      //exit(0);
          }
          break;
        case PRES:
          // --- const pressure is only defined at top well block, need to 
          // --- add gravity to other well blocks ---
          gra = 0.0;
          // for each well block
          for( i = 0; i < w[iw].wLens; i++ ) {
               ii = w[iw].wLocIJK[0][i];
               jj = w[iw].wLocIJK[1][i];
               kk = w[iw].wLocIJK[2][i];
               ijk = g.getIndex(ii, jj, kk); 
               vol = g.getBdx( ii ) * g.getBdy( jj ) * g.getBdz( kk );
               if( recalA ) diagA[ijk] -=  w[iw].wTrans[i]/vol;
               if( i > 0) 
                   gra += - fluidPtr->getGamma(ijk)
                          * (  g.getZ( kk ) 
                             - g.getZ( w[iw].wLocIJK[2][i-1] ) 
                            );
               RHS[ijk] -= w[iw].wTrans[i] * ( w[iw].wStr + gra ) / vol;
	       if(debug){
	          cout<<"Pressure Well: "
		      <<"i = "     << ii  <<' '
		      <<"j = "     << jj  <<' '
		      <<"k = "     << kk  <<' '
		      <<"ijk = "   << ijk <<' '
		      <<"diagA = " << diagA[ijk] <<' '
		      <<"RHS   = " << RHS[ijk]   <<endl;
	       }
          }
          break;
        default:
          cerr << "Well control type ERROR" << w[iw].wConType << "\n";
      }
  }
}

//==============================================================================
void P0Eqn::calcAccu(bool recalA, double dt) {
  int i;
  if( recalA ) 
      for(i = 0; i < nNode; i++) 
	  diagA[i] -= fluidPtr->getPoro_crf( i )/dt;
      
  for( i = 0; i < nNode; i++) 
       RHS[i] -=  fluidPtr->getPoro_crf( i )/dt*X[i];   
}                   

//===== ( dP0dXi and dP0dt ) =====
void P0Eqn::calcDerivatives(int iter, double dt) {
  int i, j, k, ijk, ijk1;
  if( nx > 1 ) {
      for( k = 0; k < nz; k++) {
           for( j = 0; j < ny; j++) {
                for (i = 0; i < nx - 1; i++) {
                     ijk  = g.getIndex( i, j, k );
                     ijk1 = g.getIndex(i + 1, j, k );
                     dP0dXi[0][iter][ijk] = (X[ijk1] - X[ijk]) / g.getDx(i);
                }
           }
      }
  }
  if( ny > 1) {
      for( k = 0; k < nz; k++) {
           for( i = 0; i < nx; i++) {
                for( j = 0; j < ny - 1; j++) {
                     ijk  = g.getIndex( i, j, k );        
                     ijk1 = g.getIndex( i, j + 1, k );
                     dP0dXi[1][iter][ijk] = (X[ijk1] - X[ijk]) / g.getDy(j);
                }
           }
      }
  }
  if( nz > 1) {
      for( j = 0; j < ny; j++) {
           for( i = 0; i < nx; i++) {
                for( k = 0; k < nz - 1; k++) {
                     ijk  = g.getIndex( i, j, k );
                     ijk1 = g.getIndex( i, j, k + 1 );
                     dP0dXi[2][iter][ijk] = (X[ijk1] - X[ijk]) / g.getDz(k);
                }
           }
      }
  }
  // --- calculate dP0dt ----------------------------
  if( fluidPtr->isCompressible() ) {
      for(i = 0; i < nNode; i++) dP0dt[iter][i] = (X[i]-X_old[i])/dt; 
  }
}

//===== ( OUTPUT ) =====
void P0Eqn::output( const char *fileName, int flag, const Control &c,
                    int slice_index ) {
  int i, j, k, ijk;
  double aa = (unit == FIELD)? psia_pa : 1;
  double lcoord;
  ofstream os(fileName, ios::out);
     
  os << setiosflags(ios::fixed | ios::showpoint) 
     << setprecision(5);
  // (x,y) two dimensional plot
  if(flag == 0) {
      for(k = 0; k < nz; k++) {
          for(j = 0; j < ny; j++) {
              for(i = 0; i < nx; i++) {
                  ijk = g.getIndex( i, j, k );  
                  os << g.getX(i) << "  " 
                     << g.getY(j) << "  " 
                     << X[ijk]/aa << endl;
              }
              os << endl;
          }
      }
  }             
                                       
  // (x,z) two dimensional plot
  if( flag == 1 ) {
      for(k = 0; k < nz; k++) {
          for(i = 0; i < nx; i++) {
              ijk = g.getIndex( i, 0, k ); // ijk = i + k * ny * nx;
              os << g.getX(i) << "  " << g.getZ(k) << "  " << X[ijk]/aa << endl;
          }
          os << endl;
      }       
  }        
          
  // Vertical diagonal surface plot resulting from three dimensional problem.
  // This would work provided nx = ny.
  if( flag == 2 ) {
      for( k = 0; k < nz; k++ ) {
           for(j = 0; j < ny; j++ ) {
               ijk = g.getIndex( j, j, k );
               lcoord = sqrt( g.getX(j) * g.getX(j) +
                              g.getY(j) * g.getY(j) );
               os << lcoord << "  " << g.getZ(k) << "  " << X[ijk]/aa << endl; 
           }
           os << endl;
      }                 
  }                  
                
  // Diagonal one dimensional plot resulting from two dimensional problem.
  // This would work provided nx = ny.
  if( flag == 3 ) {
      for( j = 0; j < ny; j++ ) {
           ijk = g.getIndex( j, j, 0 ); //ijk = j + j * nx;
           lcoord = sqrt( g.getX(j) * g.getX(j) +
                          g.getY(j) * g.getY(j) );
           os << lcoord << "  " << X[ijk]/aa << endl;
      }         
  }                 
                            
  // Diagonal one dimensional plot resulting from two dimensional problem.
  // This would work provided nx = nz.
  if( flag == 4 ) {
      for( k = 0; k < nz; k++ ) {
           ijk = g.getIndex( k, 0, k ); 
           lcoord = sqrt( g.getX(k) * g.getX(k) +
                          g.getZ(k) * g.getZ(k) );
           os << lcoord << "  " << X[ijk]/aa << endl;
      }                
  }                       
       
  // slice_index belongs to x-direction
  if( flag == 5 ) {
      for( k = 0; k < nz; k++ ) {
           for( j = 0; j < ny; j++ ) {
                ijk = g.getIndex( slice_index, j, k ); 
                os << g.getY(j) << "  " << g.getZ(k) << "  " << X[ijk]/aa << endl;
           }    
           os << endl;
      }
  }
                       
  // slice_index belongs to y-direction
  if( flag == 6 ) {
      for( k = 0; k < nz; k++ ) {
           for( i = 0; i < nx; i++ ) {
                ijk = g.getIndex( i, slice_index, k );
                os << g.getX(i) << "  " << g.getZ(k) << "  " << X[ijk]/aa << endl;
           }                  
           os << endl;
      }                 
  }                       
       
  // slice_index belongs to z-direction
  if( flag == 7 ) {
      for( j = 0; j < ny; j++ ) {
           for( i = 0; i < nx; i++ ) {
                ijk = g.getIndex( i, j, slice_index );
                os << g.getX(i) << "  " << g.getY(j) << "  " << X[ijk]/aa << endl;
           }           
           os << endl;
      }                
  }           
             
  // This option is used to save the solution vector in one column.
  if( flag == 8 ) {
      os << nNode << endl;
      for( i = 0; i < nNode; i++ ) os << X[i]/aa << endl;
  }

  // This option is used to save the solution vector in one column.
  if( flag == 9 ) {
      for(k = 0; k < nz; k++ ) {          
          for(j = 0; j < ny; j++ ) {
              for(i = 0; i < nx; i++ ) {
                  ijk = g.getIndex( i, j, k );
                  os << i + 1 << "  " 
                     << j + 1 << "  " 
                     << k + 1 << "  "
                     << X[ijk]/aa
                     << endl;
              }
	      os << endl;
          }
      }
  }
  os.close();
}

