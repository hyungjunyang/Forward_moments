/*
 * File: P2Eqn.cc
 * ----------------------------------
 * Implementation for P2Eqn class
 */
#include "P2Eqn.h"
#include "CYPEqn.h"
#include "P0Eqn.h"
#include "Solver.h"
#include "Control.h"
#include "Well.h"

//===== ( P2Eqn() ) =====
P2Eqn::P2Eqn( const P0Eqn &P0ee, const CYPEqn &CYPee, bool debug_ ) 
  : Eqn( P0ee ), P0e( &P0ee ), debug(debug_) {
	  
  if(debug) cout << "P2Eqn::P2Eqn() " << endl;

  int i;
  p0  = P0e->getP0(); 
  CYP = CYPee.getCYP();
  OldCYPDiag = CYPee.getOldCYPDiag();

  dP2dXi[0] = (nx > 1) ? new double[ nNode ] : NULL;
  dP2dXi[1] = (ny > 1) ? new double[ nNode ] : NULL;
  dP2dXi[2] = (nz > 1) ? new double[ nNode ] : NULL;
  for( i = 0; i < nNode; i++ ) X[i] = 0.0;
}

//===== ( ~P2Eqn() ) =====
P2Eqn::~P2Eqn() {
  if(debug) cout << "P2Eqn::~P2Eqn() " << endl;
  for( int i = 0; i < 3; i++ ) delete[] dP2dXi[i];
}

//===== ( solve() ) =====
void P2Eqn::solve( int nWells, const Well *w, const Control &c, 
                   Solver &s, int iter ) {
  // --- currently, no need to record any old value ---
  int i;
  double dt = c.getDt( iter ); bool doLU = false;
  
  // --- calculate P2 from time 0 to current time----
  // --- point to P0's derivatives for current timestep -------
  for (int it = 0; it <= iter; it++) { 
      // --- point to P0's derivatives for current timestep -------
      if( nx > 1 ) dP0dXi[0] = P0e->getDP0DXi( 0, it ); 
      if( ny > 1 ) dP0dXi[1] = P0e->getDP0DXi( 1, it ); 
      if( nz > 1 ) dP0dXi[2] = P0e->getDP0DXi( 2, it ); 
      if( fluidPtr->isCompressible() ) dP0dt = P0e->getDP0Dt( iter ); 

      // B for internal points
      calcRHS( c.getBType(), c.getBCond() );

      // boundary points
      //setBoundaries( false, c.bType, c.bCond );

      // well blocks
      setWells( false, nWells, w );

      // time-derivative
      if( fluidPtr->isCompressible() ) calcAccu( false, c.getDt( it ) );

      // --- Solve for P2 ----------------        
      s.solve( false, *this, 1, RHS );

      for( i = 0; i < nNode; i++ ) X[i] = RHS[i];

      // --- calculate P2's spatial derivatives ---------
      //calcDerivatives( c.getBType(), c.getBCond() );
  } 
}

//=====( calcRHS )==============================================================
void P2Eqn::calcRHS( const BType *bType, const double *bCond ) {
  int nxy = nx * ny;
  int i,  j,  k ;
  int i2, j2, k2;
  int i2_p1, i2_n1;
  int j2_p1, j2_n1;
  int k2_p1, k2_n1;
  int ijk, ijk_p1, ijk_n1, ijkc;
  double   cyy_p1, cyy_n1, cyy_p0;
  double   cyp_p1, cyp_n1, cyp_p0, cyp_n0;
  double perm_interface[3];
    
  for(i = 0; i < nNode; i++) RHS[i] = 0.0;        // set zero

  // --- X direction (dCYPdXidXj, i==j) -------------------
  if(nx > 1 ) {
     double dx0  = g.getDx( 0 );
     double dx02 = dx0 * dx0;
     double dxn  = g.getDx( nx - 2 );
     double dxn2 = dxn * dxn;
     double gbdx0 = g.getBdx( 0 );
     double gbdxn = g.getBdx( nx - 1 );
     for(k = 0; k < nz; k++ ) {
         //k2   = g.kPerm[k];
         k2 = g.toKPerm( k );
         for(j = 0; j < ny; j++ ) {
             //j2 = g.jPerm[j];
             j2 = g.toJPerm( j );
             for(i = 1; i < nx-1; i++) {
                 //i2     = g.iPerm[i];
                 i2     = g.toIPerm( i );
                 i2_p1  = i2 + 1;
                 i2_n1  = i2 - 1;
                 cyy_p1 = permPtr->getCYY(i2_p1, j2, k2, i2_p1, j2, k2);
                 cyy_n1 = permPtr->getCYY(i2_n1, j2, k2, i2_n1, j2, k2);
                 ijk    = g.getIndex(i  , j, k);
                 ijk_n1 = g.getIndex(i-1, j, k); 
                 permPtr->Perm_x( i, j, k, perm_interface );
                 RHS[ijk] -= 0.5*( perm_interface[2] * cyy_p1 * dP0dXi[0][ijk   ] - 
                                   perm_interface[0] * cyy_n1 * dP0dXi[0][ijk_n1]
                                 ) / g.getBdx( i );
                 ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                 ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                 cyp_p0 = CYP[ ijk_p1 * nNode + ijk  ];
                 cyp_n0 = CYP[ ijk_n1 * nNode + ijk  ];
                 ijkc   = g.getIndex(i + 1, j, k);
                 cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                 ijkc   = g.getIndex(i - 1, j, k);
                 cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                 RHS[ijk] -= ( perm_interface[2] * (cyp_p1 - cyp_p0)/g.getDx(i) 
                              -perm_interface[0] * (cyp_n0 - cyp_n1)/g.getDx(i-1)  
                             ) / g.getBdx( i );
             }
         }
     }

     if(g.getGridType() == POINT_DIS ) {
        i     = 0 ;
        //i2    = g.iPerm[ i  ];
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[0] ) {
          case CONST_PRES: 
            for(k = 0; k < nz; k++) 
               for(j = 0; j < ny; j++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ]; 
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    //j2     = g.jPerm[j];
                    j2     = g.toJPerm( j );
                    cyy_p0 = permPtr->getCYY(i2   , j2, k2, i2   , j2, k2);
                    cyy_p1 = permPtr->getCYY(i2_p1, j2, k2, i2_p1, j2, k2);
                    cyy_n1 = permPtr->getCYY(i2_p1, j2, k2, i2   , j2, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_x( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk  ];
                    ijkc   = g.getIndex(i + 1, j, k);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    RHS[ijk] -= perm_interface[2] / dx0 
                              * ( bCond[0] / perm_interface[1]
                                  *(cyy_p0+cyy_p1 -2.*cyy_n1)
                                 +cyy_p1 * dP0dXi[0][ijk]
                                 +2. * (cyp_p1 - cyp_p0) / dx0
                                );
               }
            }
            break;
        }
        i     = nx - 1 ;
        //i2    = g.iPerm[ i  ];
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[1] ) {
          case CONST_PRES: 
            for(k = 0; k < nz; k++) 
               for(j = 0; j < ny; j++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ]; 
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    //j2     = g.jPerm[j];
                    j2     = g.toJPerm( j );
                    cyy_p0 = permPtr->getCYY(i2   , j2, k2, i2   , j2, k2);
                    cyy_p1 = permPtr->getCYY(i2_n1, j2, k2, i2   , j2, k2);
                    cyy_n1 = permPtr->getCYY(i2_n1, j2, k2, i2_n1, j2, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_x( i, j, k, perm_interface );
                    ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk  ];
                    ijkc   = g.getIndex(i - 1, j, k);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc ];
                    ijk_n1 = g.getIndex(i - 1, j, k);
                    RHS[ijk] -= perm_interface[0] / dxn
                              * ( bCond[1]/perm_interface[1]*(cyy_p0+cyy_n1 -2.*cyy_p1)
                                 -cyy_n1 * dP0dXi[0][ijk_n1]
                                 +2. * (cyp_n1 - cyp_n0) / dxn
                                );
               }
            }
            break;
        }
     } else {
        i     = 0;
        //i2    = g.iPerm[ i  ];
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[0] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ];
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    //j2     = g.jPerm[j];
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i2_p1, j2, k2, i2_p1, j2, k2);
                    cyy_n1 = permPtr->getCYY(i2_p1, j2, k2, i2_n1, j2, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_x( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                    ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                    ijkc   = g.getIndex(i + 1, j, k);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    RHS[ijk] -= (0.5 * perm_interface[2] * cyy_p1 * dP0dXi[0][ijk]
                                     - perm_interface[0] * cyy_n1
                                       * (p0[ijk]-bCond[0])/gbdx0
                                ) / gbdx0;
                    RHS[ijk] -= (         perm_interface[2] * (cyp_p1 - cyp_p0) / dx0
                                 - 2.0 *  perm_interface[0] *  cyp_n0 / gbdx0
                                ) / gbdx0;
                }
            }      
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                // k2   = g.kPerm[ k  ];
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    //j2     = g.jPerm[ j  ];
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i2_p1, j2, k2, i2_p1, j2, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_x( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                    ijkc   = g.getIndex(i + 1, j, k);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    RHS[ijk] -= perm_interface[2] / gbdx0
                                *(   cyy_p1/2. * dP0dXi[0][ijk]
                                  + (cyp_p1 - cyp_p0) / dx0
                                 );
                }
            } 
            break;
        }
        i     = nx - 1;
        //i2    = g.iPerm[ i  ];
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[1] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ];
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    //j2     = g.jPerm[j];
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i2_p1, j2, k2, i2_p1, j2, k2);
                    cyy_n1 = permPtr->getCYY(i2_p1, j2, k2, i2_n1, j2, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_x( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                    ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                    ijkc   = g.getIndex(i - 1, j, k);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk_n1 = g.getIndex(i - 1, j, k);
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * (bCond[1]-p0[ijk])/gbdxn
                                 -perm_interface[0] * cyy_n1 * dP0dXi[0][ijk_n1] /2.
                                ) / gbdxn;
                    RHS[ijk] += ( 2.0 * perm_interface[2] *  cyp_p0 / gbdxn
                                      + perm_interface[0] * (cyp_n0 - cyp_n1) / dxn
                                  ) / gbdxn;
                }
            }      
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ];
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    //j2     = g.jPerm[j];
                    j2     = g.toJPerm( j );
                    cyy_n1 = permPtr->getCYY(i2_n1, j2, k2, i2_n1, j2, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_x( i, j, k, perm_interface );
                    ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                    ijkc   = g.getIndex(i - 1, j, k);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk_n1 = g.getIndex(i - 1, j, k);    
                    RHS[ijk] -=  perm_interface[0]/gbdxn 
                                *(-cyy_n1 /2. * dP0dXi[0][ijk_n1]
                                  -(cyp_n0 - cyp_n1) / dxn
                                 );
                }
            } 
            break;
        }
     }
  }
  // --- Y direction (dCYPdXidXj, i==j) -------------------
  if(ny > 1) {
     double dy0   = g.getDy( 0 );
     double dy02  = dy0 * dy0;
     double dyn   = g.getDy( ny - 2 );
     double gbdy0 = g.getBdy( 0 );
     double gbdyn = g.getBdy( ny - 1 );
     for(k = 0; k < nz; k++) {
         //k2   = g.kPerm[k];
         k2 = g.toKPerm( k );
         for(i = 0; i < nx; i++) {
             //i2   = g.iPerm[i];
             i2   = g.toIPerm( i );
             for(j = 1; j < ny - 1; j++) {
                 //j2     = g.jPerm[j];
                 j2     = g.toJPerm( j );
                 j2_p1  = j2 + 1;
                 j2_n1  = j2 - 1;
                 cyy_p1 = permPtr->getCYY(i2, j2_p1, k2, i2, j2_p1, k2);
                 cyy_n1 = permPtr->getCYY(i2, j2_n1, k2, i2, j2_n1, k2);
                 ijk    = g.getIndex(i, j  , k);
                 ijk_n1 = g.getIndex(i, j-1, k); 
                 permPtr->Perm_y( i, j, k, perm_interface );
                 RHS[ijk] -= 0.5*( perm_interface[2] * cyy_p1 * dP0dXi[1][ijk   ] - 
                                   perm_interface[0] * cyy_n1 * dP0dXi[1][ijk_n1]
                                 ) / g.getBdy( j );
                 ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                 ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                 cyp_p0 = CYP[ ijk_p1 * nNode + ijk  ];
                 cyp_n0 = CYP[ ijk_n1 * nNode + ijk  ];
                 ijkc   = g.getIndex(i, j + 1, k);
                 cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                 ijkc   = g.getIndex(i, j - 1, k);
                 cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                 RHS[ijk] -= ( perm_interface[2] * (cyp_p1 - cyp_p0)/g.getDy( j ) 
                              -perm_interface[0] * (cyp_n0 - cyp_n1)/g.getDy(j-1)  
                             ) / g.getBdy( j );
             }
         }
     }
     if(g.getGridType() == POINT_DIS ) {
        j     = 0 ;
        //j2    = g.jPerm[ j ];
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        j2_n1 = j2 - 1;
        switch( bType[2] ) {
          case CONST_PRES: 
            for(k = 0; k < nz; k++) 
               for( i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ]; 
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    //i2     = g.iPerm[ i  ];
                    i2     = g.toIPerm( i );
                    cyy_p0 = permPtr->getCYY(i2, j2   , k2, i2, j2   , k2);
                    cyy_p1 = permPtr->getCYY(i2, j2_p1, k2, i2, j2_p1, k2);
                    cyy_n1 = permPtr->getCYY(i2, j2_p1, k2, i2, j2   , k2);
                    ijk    = g.getIndex(i, j, k);
                    ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk  ];
                    ijkc   = g.getIndex(i, j + 1, k);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] / dy0 
                              * ( bCond[2] / perm_interface[1]*(cyy_p0+cyy_p1 -2.*cyy_n1)
                                 +cyy_p1 * dP0dXi[1][ijk]
                                 +2. * (cyp_p1 - cyp_p0) / dy0
                                );
               }
            }
            break;
        }
        j     = ny - 1 ;
        //j2    = g.jPerm[ j  ];
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        j2_n1 = j2 - 1;
        switch( bType[3] ) {
          case CONST_PRES: 
            for(k = 0; k < nz; k++) 
               for(i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ]; 
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    //i2     = g.iPerm[ i  ];
                    i2     = g.toIPerm( i );
                    cyy_p0 = permPtr->getCYY(i2, j2   , k2, i2, j2   , k2);
                    cyy_p1 = permPtr->getCYY(i2, j2_n1, k2, i2, j2   , k2);
                    cyy_n1 = permPtr->getCYY(i2, j2_n1, k2, i2, j2_n1, k2);
                    ijk    = g.getIndex(i, j, k);
                    ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk  ];
                    ijkc   = g.getIndex(i, j - 1, k);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc ];
                    ijk_n1 = g.getIndex(i, j - 1, k);
                    permPtr->Perm_y( i, j, k, perm_interface ); 
                    RHS[ijk] -= perm_interface[0] / dyn
                              * ( bCond[3]/perm_interface[1]*(cyy_p0+cyy_n1 -2.*cyy_p1)
                                 -cyy_n1 * dP0dXi[1][ijk_n1]
                                 +2. * (cyp_n1 - cyp_n0) / dyn
                                );
               }
            }
            break;
        }
     } else {
        j     = 0 ;
        //j2    = g.jPerm[ j  ];
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        j2_n1 = j2 - 1;
        switch( bType[2] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ];
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    //i2     = g.iPerm[ i  ];
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i2, j2_p1, k2, i2, j2_p1, k2);
                    cyy_n1 = permPtr->getCYY(i2, j2_n1, k2, i2, j2_n1, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                    ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                    ijkc   = g.getIndex(i, j + 1, k);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    RHS[ijk] -= (0.5 * perm_interface[2] * cyy_p1 * dP0dXi[1][ijk]
                                     - perm_interface[0] * cyy_n1 * (p0[ijk]-bCond[2])/gbdy0
                                ) / gbdy0;
                    RHS[ijk] -= (         perm_interface[2] * (cyp_p1 - cyp_p0) / dy0
                                 - 2.0 *  perm_interface[0] *  cyp_n0 / gbdy0
                                ) / gbdy0;
                }
            }      
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ];
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    //i2     = g.iPerm[ i  ];
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i2, j2_p1, k2, i2, j2_p1, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                    ijkc   = g.getIndex(i, j + 1, k);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    RHS[ijk] -= perm_interface[2] / gbdy0
                                *(   cyy_p1/2. * dP0dXi[1][ijk]
                                  + (cyp_p1 - cyp_p0) / dy0
                                 );
                }
            }
            break;
        }
        j     = ny - 1;
        //j2    = g.jPerm[ j  ];
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        j2_n1 = j2 - 1;
        switch( bType[3] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ];
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    //i2     = g.iPerm[ i  ];
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i2, j2_p1, k2, i2, j2_p1, k2);
                    cyy_n1 = permPtr->getCYY(i2, j2_n1, k2, i2, j2_n1, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                    ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                    ijkc   = g.getIndex(i, j - 1, k);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk_n1 = g.getIndex(i, j - 1, k);
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * (bCond[3]-p0[ijk])/gbdyn
                                 -perm_interface[0] * cyy_n1 * dP0dXi[1][ijk_n1] /2.
                                ) / gbdyn;
                    RHS[ijk] += ( 2.0 * perm_interface[2] *  cyp_p0 / gbdyn
                                      + perm_interface[0] * (cyp_n0 - cyp_n1) / dyn
                                  ) / gbdyn;
                }
            }
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                //k2   = g.kPerm[ k  ];
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    //i2     = g.iPerm[ i  ];
                    i2     = g.toIPerm( i );
                    cyy_n1 = permPtr->getCYY(i2, j2_n1, k2, i2, j2_n1, k2);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                    ijkc   = g.getIndex(i, j - 1, k);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk_n1 = g.getIndex(i, j - 1, k);    
                    RHS[ijk] -=  perm_interface[0]/gbdyn 
                                *(-cyy_n1 /2. * dP0dXi[1][ijk_n1]
                                  -(cyp_n0 - cyp_n1) / dyn
                                 );
                }
            }
            break;
        }
     }
  }
  // --- Z direction (dCYPdXidXj, i==j) -------------------
  if(nz > 1) {
     double dz0   = g.getDz( 0 );
     double dz02  = dz0 * dz0;
     double dzn   = g.getDz( nz - 2 );
     double gbdz0 = g.getBdz( 0 );
     double gbdzn = g.getBdz( nz - 1 );
     for(j = 0; j < ny; j++) {
         //j2 = g.jPerm[ j  ];
         j2 = g.toJPerm( j );
         for(i = 0; i < nx; i++) {
             //i2 = g.iPerm[i];
             i2 = g.toIPerm( i );
             for(k = 1; k < nz - 1; k++) {
                 //k2   = g.kPerm[k];
                 k2     = g.toKPerm( k );
                 k2_p1  = k2 + 1;
                 k2_n1  = k2 - 1;
                 cyy_p1 = permPtr->getCYY(i2, j2, k2_p1, i2, j2, k2_p1);
                 cyy_n1 = permPtr->getCYY(i2, j2, k2_n1, i2, j2, k2_n1);
                 ijk    = g.getIndex(i, j, k    );
                 ijk_n1 = g.getIndex(i, j, k - 1); 
                 permPtr->Perm_z( i, j, k, perm_interface );
                 RHS[ijk] -= 0.5*( perm_interface[2] * cyy_p1 * dP0dXi[2][ijk   ] - 
                                   perm_interface[0] * cyy_n1 * dP0dXi[2][ijk_n1]
                                 ) / g.getBdz( k );
                 ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                 ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                 cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                 cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                 ijkc   = g.getIndex(i, j, k + 1);
                 cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                 ijkc   = g.getIndex(i, j, k - 1);
                 cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                 RHS[ijk] -= ( perm_interface[2] * (cyp_p1 - cyp_p0)/g.getDz(k) 
                              -perm_interface[0] * (cyp_n0 - cyp_n1)/g.getDz(k-1)  
                             ) / g.getBdz( k );
                 RHS[ijk] -= 0.5 * fluidPtr->getGamma(ijk) 
                                 * ( perm_interface[2] * cyy_p1 
                                    -perm_interface[0] * cyy_n1 ) / g.getBdz( k);
             }
         }
     }
     if(g.getGridType() == POINT_DIS ) {
        k     = 0 ;
        //k2    = g.kPerm[ k  ];
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[4] ) {
          case CONST_PRES:
            for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k) ] = 0.0;
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++ ) {
                //j2 = g.jPerm[j];
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    //i2   = g.iPerm[ i  ]; 
                    i2     = g.toIPerm( i );
                    cyy_p0 = permPtr->getCYY(i2, j2, k2   , i2, j2, k2   );
                    cyy_p1 = permPtr->getCYY(i2, j2, k2_p1, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i2, j2, k2_p1, i2, j2, k2   );
                    ijk    = g.getIndex(i, j, k);
                    ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk  ];
                    ijkc   = g.getIndex(i, j, k + 1);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] / dz0 
                              * ( bCond[4] /perm_interface[1]*(cyy_p0+cyy_p1 -2.*cyy_n1)
                                 +cyy_p1 * ( dP0dXi[2][ijk] + fluidPtr->getGamma(ijk) )
                                 +2. * (cyp_p1 - cyp_p0) / dz0
                                );
                }
            }
            break;
        }
        k     = nz - 1 ;
        //k2    = g.kPerm[ k  ];
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[5] ) {
          case CONST_PRES: 
            for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k) ] = 0.0;
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++ ) {
                //j2   = g.jPerm[j];
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    //i2   = g.iPerm[ i  ]; 
                    i2     = g.toIPerm( i );
                    cyy_p0 = permPtr->getCYY(i2, j2, k2   , i2, j2, k2   );
                    cyy_p1 = permPtr->getCYY(i2, j2, k2_n1, i2, j2, k2   );
                    cyy_n1 = permPtr->getCYY(i2, j2, k2_n1, i2, j2, k2_n1);
                    ijk    = g.getIndex(i, j, k);
                    ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk  ];
                    ijkc   = g.getIndex(i, j, k - 1);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc ];
                    ijk_n1 = g.getIndex(i, j, k - 1);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] / dzn
                              * ( bCond[5]/perm_interface[1]*(cyy_p0+cyy_n1 -2.*cyy_p1)
                                 -cyy_n1 * ( dP0dXi[2][ijk_n1] + fluidPtr->getGamma(ijk) )
                                 +2. * (cyp_n1 - cyp_n0) / dzn
                                );
               }
            }
            break;
        }
     } else {
        k     = 0 ;
        //k2    = g.kPerm[ k  ];
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[4] ) {
          case CONST_PRES:
            for(j = 0; j < ny; j++ ) {
                //j2   = g.jPerm[j];
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    //i2   = g.iPerm[ i  ]; 
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i2, j2, k2_p1, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i2, j2, k2_n1, i2, j2, k2_n1);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                    ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                    ijkc   = g.getIndex(i, j, k + 1);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    RHS[ijk] -= (0.5 * perm_interface[2] * cyy_p1 
                                     * (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                     - perm_interface[0] * cyy_n1 
                                     * (  (p0[ijk]-bCond[4])/gbdz0 
                                        + fluidPtr->getGamma(ijk)/2.)
                                ) / gbdz0;
                    RHS[ijk] -= (         perm_interface[2] * (cyp_p1 - cyp_p0) / dz0
                                 - 2.0 *  perm_interface[0] *  cyp_n0 / gbdz0
                                ) / gbdz0;
                }
            }
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++ ) {
                //j2   = g.jPerm[j];
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    //i2   = g.iPerm[ i  ]; 
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i2, j2, k2_p1, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i2, j2, k2_n1, i2, j2, k2_n1);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                    ijkc   = g.getIndex(i, j, k + 1);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    RHS[ijk] -= perm_interface[2] / gbdz0
                                *(  cyy_p1/2. * ( dP0dXi[2][ijk]
                                                 +fluidPtr->getGamma(ijk)
                                                )
                                  + (cyp_p1 - cyp_p0) / dz0
                                  +  cyy_n1/2. * fluidPtr->getGamma(ijk)
                                 );
                    RHS[ijk] +=  perm_interface[0] / gbdz0 
                                *cyy_n1/2 * fluidPtr->getGamma( ijk );
                }
            }
            break;
        }
        k     = nz - 1;
        //k2    = g.kPerm[ k  ];
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[5] ) {
          case CONST_PRES:
            for(j = 0; j < ny; j++ ) {
                //j2   = g.jPerm[j];
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    //i2   = g.iPerm[ i  ]; 
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i2, j2, k2_p1, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i2, j2, k2_n1, i2, j2, k2_n1);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                    ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                    cyp_p0 = CYP[ ijk_p1 * nNode + ijk ];
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                    ijkc   = g.getIndex(i, j, k - 1);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk_n1 = g.getIndex(i, j, k - 1);
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 
                                  *( (bCond[5]-p0[ijk])/gbdzn + fluidPtr->getGamma(ijk)/2.)
                                 -perm_interface[0] * cyy_n1 
                                  *( dP0dXi[2][ijk_n1] + fluidPtr->getGamma(ijk)) /2.
                                ) / gbdzn;
                    RHS[ijk] += ( 2.0 * perm_interface[2] *  cyp_p0 / gbdzn
                                      + perm_interface[0] * (cyp_n0 - cyp_n1) / dzn
                                ) / gbdzn;
                }
            }
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++ ) {
                //j2   = g.jPerm[j];
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    //i2   = g.iPerm[ i  ]; 
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i2, j2, k2_p1, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i2, j2, k2_n1, i2, j2, k2_n1);
                    ijk    = g.getIndex(i, j, k);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                    cyp_n0 = CYP[ ijk_n1 * nNode + ijk ];
                    ijkc   = g.getIndex(i, j, k - 1);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk_n1 = g.getIndex(i, j, k - 1);    
                    RHS[ijk] -=  perm_interface[0]/gbdzn 
                                *(-cyy_n1 /2. * (dP0dXi[2][ijk_n1]+fluidPtr->getGamma(ijk))
                                  -cyy_p1 /2. * fluidPtr->getGamma(ijk)
                                  -(cyp_n0 - cyp_n1) / dzn
                                 );
                    RHS[ijk] -=  perm_interface[2]/gbdzn*cyy_p1 /2.
                               * fluidPtr->getGamma(ijk);
                }
            }
            break;
        }
     }
  }
}

//============================================================================
void P2Eqn::setWells(bool recalA, int nWells, const Well *w)
{
  int i, ijk, iw, ii, jj, kk;
  double vol, gra;
  for(iw = 0; iw < nWells; iw++ ) {
      switch( w[iw].wConType ) {
        // constant Rate
        case RATE:
          // --- only deal with single block penetration now (AAAA) ---
          // for each well block
          for( i = 0; i < w[iw].wLens; i++ ) {
              //ijk,vol: well block index and volume
              ii = w[iw].wLocIJK[0][i];
              jj = w[iw].wLocIJK[1][i];
              kk = w[iw].wLocIJK[2][i];
              ijk = ii + ( jj+kk*ny )*nx;
              vol = g.getBdx( ii ) * g.getBdy( jj ) * g.getBdz(kk);
              RHS[ijk] -= 0.5 * regnPtr->getNodePermVar(ijk) * 
                (w[iw].wStr/vol);
          }
          break;
        // contant pressure
        case PRES:
          // --- const pressure is only dePermd at top well block, need to 
          // --- add gravity to other well blocks ---
          gra = 0.0;
          // for each well block
          for(i = 0; i < w[iw].wLens; i++) {
              ii = w[iw].wLocIJK[0][i];
              jj = w[iw].wLocIJK[1][i];
              kk = w[iw].wLocIJK[2][i];              
              ijk = ii + ( jj+kk*ny ) * nx;      
              vol = g.getBdx( ii) * g.getBdy(jj) * g.getBdz(kk);
              RHS[ijk] -= w[iw].wTrans[i] * CYP[ijk*nNode+ijk] / 
                vol;
              if(i > 0) gra += - fluidPtr->getGamma(ijk)
                        * ( g.getZ( kk ) - g.getZ( w[iw].wLocIJK[2][i-1] ) 
                          );
              RHS[ijk] -= 0.5 * regnPtr->getNodePermVar(ijk) *
                w[iw].wTrans[i] * ( -p0[ijk]+(w[iw].wStr+gra) ) / vol;
          }
          break;
        default:
          cerr << "Well control type ERROR" << w[iw].wConType << "\n";
      }
  }
}

//===== ( calcAccu ) =====
void P2Eqn::calcAccu(bool recalA, double dt) {
  for(int k=0; k<nz; k++) {
      for(int j=0; j<ny; j++) {
          for(int i=0; i<nx; i++) {
              int ijk = i+(j+k*ny)*nx;
              //  before solve, X(P2) still store the old timestep value
              RHS[ijk] -=  fluidPtr->getPoro_crf( ijk ) / dt * X[ijk];
              RHS[ijk] -=  fluidPtr->getPoro_crf( ijk )
                          *( CYP[ijk*(1+nNode)] - OldCYPDiag[ijk]
                            )/dt;
              RHS[ijk] +=  0.5* fluidPtr->getPoro_crf( ijk )
                          * regnPtr->getNodePermVar(ijk) * dP0dt[ijk];
          }
      }
  }
}

//============================================================================
void P2Eqn::calcDerivatives( const BType *bType, const double *bCond  )   
{
  int i, j, k, ijk, nxy;
  cout << "P2Eqn::calcDerivatives() is not available!"<<endl;
/*
  // --- calculate dP2dXi[0] ------------------- 
  if(nx>1) {
    for(k=0; k<nz; k++) for(j=0; j<ny; j++) {
      ijk = (j+k*ny)*nx;
      dP2dXi[0][ijk] = (g.gridType==POINT_DIS&&bType[0]==CONST_RATE)? 
        -0.5* permPtr->CY[ijk+ijk*nNode]*bCond[0]/perm_interface[1] :
        (X[ijk+1]-X[ijk])/g.dx[0];
     
      for (i=1; i<nx-1; i++) {
        ijk=i+(j+k*ny)*nx;
        dP2dXi[0][ijk] = (X[ijk+1]-X[ijk-1])/(g.dx[i]+g.dx[i-1]);
      }

      ijk += 1;
      dP2dXi[0][ijk] = (g.gridType==POINT_DIS&&bType[1]==CONST_RATE)? 
        0.5* permPtr->CY[ijk+ijk*nNode]*bCond[1]/perm_interface[1] :
        (X[ijk]-X[ijk-1])/g.dx[nx-2];
    }
  }
  
  // --- calculate dP2dXi[1] ------------------- 
  if(ny>1) {
    for(k=0; k<nz; k++) for(i=0; i<nx; i++) {
      ijk = i+k*ny*nx;
      dP2dXi[1][ijk] = (g.gridType==POINT_DIS&&bType[2]==CONST_RATE) ? 
        -0.5* permPtr->CY[ijk+ijk*nNode]*bCond[2]/perm_interface[1] :
        (X[ijk+nx]-X[ijk])/g.dy[0];
    
      for (j=1; j<ny-1; j++) {
        ijk=i+(j+k*ny)*nx;
        dP2dXi[1][ijk] = (X[ijk+nx]-X[ijk-nx])/(g.dy[j]+g.dy[j-1]);
      } 

      ijk += nx;
      dP2dXi[1][ijk] = (g.gridType==POINT_DIS&&bType[3]==CONST_RATE)? 
        0.5* permPtr->CY[ijk+ijk*nNode]*bCond[3]/perm_interface[1] :
        (X[ijk]-X[ijk-nx])/g.dy[ny-2];
    }
  }

  // --- calculate dP2dXi[2] ------------------- 
  if(nz>1) {
    nxy = nx*ny;
    for(j=0; j<ny; j++) for(i=0; i<nx; i++) {
      ijk = i+j*nx;
      dP2dXi[2][ijk] = (g.gridType==POINT_DIS&&bType[4]==CONST_RATE)? 
        -0.5* permPtr->CY[ijk+ijk*nNode]*bCond[4]/perm_interface[1] :
        (X[ijk+nxy]-X[ijk])/g.dz[0];
    
      for (k=1; k<nz-1; k++) {
        ijk=i+(j+k*ny)*nx;
        dP2dXi[2][ijk] = (X[ijk+nxy]-X[ijk-nxy])/(g.dz[k]+g.dz[k-1]);
      } 

      ijk += nxy;
      dP2dXi[2][ijk] = (g.gridType==POINT_DIS&&bType[5]==CONST_RATE)? 
        0.5* permPtr->CY[ijk+ijk*nNode]*bCond[5]/perm_interface[1] :
        (X[ijk]-X[ijk-nxy])/g.dz[nz-2];
    }
  }
  */
}

//===== ( Output() ) =====
void P2Eqn::output( const char *fileName, int flag, const Control &c,
                    int slice_index ) { 
  int i, j, k, ijk;
  double aa = (unit==FIELD)? psia_pa : 1;
  ofstream os(fileName, ios::out);
  double lcoord;

  os << setiosflags(ios::fixed | ios::showpoint)
     << setprecision(5);

  // (x,y) two dimensional plot
  if( flag == 0 ) {
      for(k = 0; k < nz; k++) {
          for(j = 0; j < ny; j++) {
              for(i = 0; i < nx; i++) {
                  ijk = i+(j+k*ny)*nx;
                  os << g.getX( i ) << "  " 
                     << g.getY( j ) << "  "
                     << ( p0[ijk]+X[ijk] ) / aa << endl;
              }
              os << endl;
          }
      }
  }
  
  // (x,z) two dimensional plot
  if(flag == 1 ) {
      for(k = 0; k < nz; k++) {
          for(i = 0; i < nx; i++) {
              ijk = i+k*ny*nx;
              os << g.getX( i ) << "  " 
                 << g.getZ( k ) << "  "
                 << ( p0[ijk]+X[ijk] ) / aa << endl;
          }
          os << endl;
      }
  }

  // Vertical diagonal surface plot resulting from three dimensional problem.
  // This would work provided nx=ny.
  if(flag == 2 ) {
      for(k = 0; k < nz; k++ ) {
          for(j = 0; j < ny; j++ ) {
              ijk = j+(j+k*ny)*nx;              
              lcoord = sqrt(  g.getX( j ) * g.getX( j ) 
                            + g.getY( j ) * g.getY( j ) );
              os << lcoord      << "  " 
                 << g.getZ( k ) << "  "
                 << ( p0[ijk]+X[ijk] ) / aa << endl; 
          }
          os << endl;
      }
  }

  // Diagonal one dimensional plot resulting from two dimensional problem.
  // This would work provided nx=ny.
  if(flag == 3 ) {
      for(j = 0; j < ny; j++ ) {
          ijk = j+j*nx;
          lcoord = sqrt(  g.getX( j ) * g.getX( j ) 
                        + g.getY( j ) * g.getY( j ) );
          os << lcoord << "  " 
             << ( p0[ijk]+X[ijk] ) / aa << endl;
      }
  }

  // Diagonal one dimensional plot resulting from two dimensional problem.
  // This would work provided nx=nz.
  if(flag == 4 ) {
      for(k = 0; k < nz; k++ ) {
          ijk = k+k*nx;
          lcoord = sqrt(  g.getX( k ) * g.getX( k ) 
                        + g.getZ( k ) * g.getZ( k ) );
          os << lcoord << "  " 
             << ( p0[ijk]+X[ijk] ) / aa << endl;
      }
  }

  // slice_index belongs to x-direction
  if(flag == 5 ) {
      for(k = 0; k < nz; k++ ) {
          for(j = 0; j < ny; j++ ) {
              ijk = slice_index + (j+k*ny)*nx;
              os << g.getY( j ) << "  " 
                 << g.getZ( k ) << "  "
                 << ( p0[ijk]+X[ijk] ) / aa << endl;
          }
          os << endl;
      }
  }

  // slice_index belongs to y-direction
  if(flag == 6 ) {
      for(k = 0; k < nz; k++ ) {
          for(i = 0; i < nx; i++ ) {
              ijk = i + (slice_index+k*ny)*nx;
              os << g.getX( i ) << "  " 
                 << g.getZ( k ) << "  "
                 << ( p0[ijk]+X[ijk] ) / aa << endl;
          }
          os << endl;
      }
  }

  // slice_index belongs to z-direction
  if(flag == 7 ) {
      for(j = 0; j < ny; j++ ) {
          for(i = 0; i < nx; i++ ) {
              ijk = i + (j+slice_index*ny)*nx;
              os << g.getX( i ) << "  " 
                 << g.getY( j ) << "  "
                 << ( p0[ijk]+X[ijk] ) / aa << endl;
          }
          os << endl;
      }
  }

  // This option is used to save the solution vector in one column.
  if(flag == 8 ) {
      os << nNode << endl;
      for(i = 0; i < nNode; i++ ) 
        os << ( p0[i]+X[i] ) / aa << endl;
  }

  os.close();
}
