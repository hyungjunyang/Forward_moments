/*
 * File: CYPEqn.cc
 * ----------------------------------
 * Implementation for CYPEqn class
 */
#include "CYPEqn.h"
#include "P0Eqn.h"
#include "Well.h"
#include "Solver.h"
#include "Control.h"

//===== CYPEqn( const P0Eqn &P0ee ) =====
CYPEqn::CYPEqn( const P0Eqn &P0ee, bool debug_ ) 
: Eqn( P0ee ), P0e( &P0ee ), debug(debug_) {
  
  if(debug) cout << "CYPEqn::CYPEqn( )"<<endl;
  
  p0  = P0e->getP0();
  
  if(debug) cout << nNodePerm <<' '<< nNode << endl;
  
  try {
    CYP = new double[ nNodePerm * nNode];
  }
  catch (bad_alloc){
    cerr << "problem in creating CYP..\n";
    exit(0);
  }
  
  for(long i = 0; i < nNodePerm * nNode; i++) {
      CYP[i] = 0.0;
  }
}

//===== ~CYPEqn() =====
CYPEqn::~CYPEqn() {
  if(debug) cout << "CYPEqn::~CYPEqn( )"<<endl;
  delete[] CYP;
}

//===== solve() =====
void CYPEqn::solve( int nWells, const Well *w, const Control &c, 
                    Solver &s, int iter ) {
       
  // --- point to P0's derivatives for current timestep -------
  if( nx > 1 ) dP0dXi[0] = P0e->getDP0DXi( 0, iter ) ;
  if( ny > 1 ) dP0dXi[1] = P0e->getDP0DXi( 1, iter ) ;
  if( nz > 1 ) dP0dXi[2] = P0e->getDP0DXi( 2, iter ) ;
  if( fluidPtr->isCompressible() ) dP0dt = P0e->getDP0Dt( iter ); 
  // --- assign X_old (if necessary) (it's the CYP(i,i)) ---
  // --- used to calculate dCYP(i,i)/dt --------------------
  // if inCompressible, then no time dependent
  if(fluidPtr->isCompressible() ) 
     for(int i = 0; i < nNode; i++ ) 
        X_old[i] = CYP[i + i * nNode];
  
  for(int kc = 0; kc < nzPerm; kc++) {
      for(int jc = 0; jc < nyPerm; jc++) {
          for(int ic = 0; ic < nxPerm; ic++) { 
              int ijkc = g.getPermIndex(ic, jc, kc);
              calcRHS(ic, jc, kc, c.getBType(), c.getBCond());
              if(fluidPtr->isCompressible() ) 
                 calcAccu( false, c.getDt( iter ) );
             
              for(int i = 0; i < nNode; i++) {
                  CYP[i + nNode * ijkc  ] = RHS[i]; 
              }
          }
      }
  }
          
  // --- Solve for CYP ----------------
  // --- calculate RHS for each reference point, combine them --- 
  // --- and solve them together by a multi-rhs solver ---------
  // s.solve( false, *this, nA, CYP );   // X and B share space in CYP
  //for(int i = 0; i < nNode; i++) {
  //    cout <<" CYP_RHS[ijkc = 0] = "  << CYP[i + 0 * nNode ] << endl;
  //}

  s.solve( false, *this, nNodePerm, CYP );

  //checkBalance();
  /*
  ofstream os("CYP_1p.out", ios::out);
  for(int i = 0; i < nNodePerm * nNode; ++i) {
      os << CYP[i] <<endl;
  } 
  os.close();
  */

  if(debug) printCYP();
}

//===== ( calcRHS() ) =====
void CYPEqn::calcRHS( int & i1_perm, int & j1_perm, int & k1_perm, 
                      const BType *bType, const double *bCond) {
  int nxy = nx * ny;
  int i , j , k ;
  int i2, j2, k2;
  int ijk, ijk_p1, ijk_n1;
  int i2_p1, i2_n1;
  int j2_p1, j2_n1;
  int k2_p1, k2_n1;

  double   cyy_p1, cyy_n1, cyy_p0;
  double perm_interface[3];
  
  for( i = 0; i < nNode; i++ ) RHS[i] = 0.0;
  
  // --- X direction (dCYPdx) -------------------------
  if(nx > 1 ) {
     double dx0   = g.getDx(0);
     double dx02  = dx0 * dx0;
     double dxn   = g.getDx( nx - 2 );
     double gbdx0 = g.getBdx( 0 );
     double gbdxn = g.getBdx( nx - 1);
     for(k = 0; k < nz; k++ ) {
         k2 = g.toKPerm( k );
         for(j = 0; j < ny; j++ ) {
             j2 = g.toJPerm( j );
             for(i = 1; i < nx - 1; i++ ) {
                 i2     = g.toIPerm( i );
                 i2_p1  = i2 + 1;
                 i2_n1  = i2 - 1;
                 ijk_n1 = g.getIndex( i-1, j, k );
                 ijk    = g.getIndex( i  , j, k );
                 ijk_p1 = g.getIndex( i+1, j, k );
                 cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                 cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                 permPtr->Perm_x( i, j, k, perm_interface );
                 RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[0][ijk   ] -
                               perm_interface[0] * cyy_n1 * dP0dXi[0][ijk_n1]
                             ) / g.getBdx( i ); 
             }
         }
     }
     if(g.getGridType() == POINT_DIS ) {
        i      = 0 ;
        i2     = g.toIPerm( i );
        i2_p1   = i2 + 1;
        i2_n1   = i2 - 1;
        switch( bType[0] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++)
                for(j = 0; j < ny; j++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2 = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2   , j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= 2.*( cyy_p1 * dP0dXi[0][ijk]
                                   + bCond[0] * (cyy_p1 - cyy_p0)/perm_interface[1]
                                   ) / dx0 * perm_interface[2];
                }
            }
            break;
        }
        i      = nx - 1 ;
        i2     = g.toIPerm( i );
        i2_p1  = i2 + 1;
        i2_n1  = i2 - 1;
        switch( bType[1] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++) 
                for(j = 0; j < ny; j++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++) {
                    j2 = g.toJPerm( j );
                    ijk    = g.getIndex(     i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2   , j2, k2);
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= 2.* perm_interface[0] / dxn 
                                  * ( cyy_n1 * ( -dP0dXi[0][ijk_n1] )
                                    - bCond[1] * (cyy_p0 - cyy_n1 )/ perm_interface[1]
                                    );
                }
            }
            break;
        }
      }
     else { // cell-centered
        i      = 0 ;
        i2     = g.toIPerm( i );
        i2_p1  = i2 + 1;
        i2_n1  = i2 - 1;
        switch( bType[0] ) {
          case CONST_PRES: //ok
            //cout<< "CONST_PRES: i = 0 " << endl;		  
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[0][ijk]
                                 -perm_interface[0] * cyy_n1 *
                                    (p0[ijk] - bCond[0])/(0.5*gbdx0)
                                ) / gbdx0;
                }
            }
            break;
          case CONST_RATE: //ok
	    //cout<< "CONST_RATE: i = 0 " << endl;
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * cyy_p1 / gbdx0 * dP0dXi[0][ijk];
                }
            }
            break;
        }

        i     = nx - 1;
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[1] ) {
          case CONST_PRES:  //ok
            //cout<< "CONST_PRES: i = nx - 1 " << endl;		  
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 *
                                   (bCond[1] - p0[ijk])/(0.5 * gbdxn)
                                 -perm_interface[0] * cyy_n1 * ( dP0dXi[0][ijk_n1])
                                ) / gbdxn;
                }
            }
            break;
          case CONST_RATE: //ok
	    //cout<< "CONST_RATE: i = nx - 1 " << endl;
	    //cout << i << ' ' << i2 << ' '<< gbdxn << endl;
	    for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    ijk    = g.getIndex( i    , j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
		    //cout << j <<' '<<dP0dXi[0][ijk_n1]<< endl;
                    RHS[ijk] -= perm_interface[0] * cyy_n1 / gbdxn * (-dP0dXi[0][ijk_n1]);
                }
            }
            break;
        }
     }
     //exit(0);
  }

  if(ny > 1 ) {
     double dy0   = g.getDy( 0 );
     double dy02  = dy0 * dy0;
     double dyn   = g.getDy( ny - 2 );
     double gbdy0 = g.getBdy( 0 );
     double gbdyn = g.getBdy( ny - 1 );
     for(k = 0; k < nz; k++ ) {
         k2 = g.toKPerm( k );
         for(i = 0; i < nx; i++ ) {
             i2 = g.toIPerm( i );
             for(j = 1; j < ny - 1; j++ ) {
                 j2     = g.toJPerm( j );
                 j2_p1  = j2 + 1;
                 j2_n1  = j2 - 1;
                 ijk_n1 = g.getIndex( i, j-1, k );
                 ijk    = g.getIndex( i, j  , k );
                 ijk_p1 = g.getIndex( i, j+1, k );
                 cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                 cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                 permPtr->Perm_y( i, j, k, perm_interface );
                 RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[1][ijk    ]
                              - 
                               perm_interface[0] * cyy_n1 * dP0dXi[1][ijk_n1 ]
                             ) / g.getBdy( j );
             }
         }
     }
     if(g.getGridType() == POINT_DIS ) {
        j      = 0 ;
        j2     = g.toJPerm( j );
        j2_p1  = j2 + 1;
        switch( bType[2] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++)
                for(i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j  , k );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2   , k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[2] /dy0 
                                   * ( cyy_p1 * dP0dXi[1][ijk]
                                      +bCond[2]* (cyy_p1 - cyy_p0 )/perm_interface[1]
                                     );
                }
            }
            break;
        }
        j      = ny - 1 ;
        j2     = g.toJPerm( j );
        j2_n1  = j2 - 1;
        switch( bType[3] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++)
                for( i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2   , k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[0] / dyn
                                   * ( cyy_n1 * ( -dP0dXi[1][ijk_n1] )
                                     - bCond[3]* (cyy_p0 - cyy_n1)/perm_interface[1]
                                     );
                }
            }
            break;
        }
     } 
     else {
        j      = 0 ;
        j2     = g.toJPerm( j );
        j2_p1  = j2 + 1;
        j2_n1  = j2 - 1;
        switch( bType[2] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[1][ijk]
                                 -perm_interface[0] * cyy_n1 *
                                    (p0[ijk] - bCond[2])/(0.5*gbdy0)
                                ) / gbdy0;
                }
            }
            break;
          case CONST_RATE: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * cyy_p1 / gbdy0 * dP0dXi[1][ijk];
                }
            }
            break;
        }
        j      = ny - 1;
        j2     = g.toJPerm( j );
        j2_p1  = j2 + 1;
        j2_n1  = j2 - 1;
        switch( bType[3] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 *
                                    (bCond[3] - p0[ijk])/(0.5*gbdyn)
                                 -perm_interface[0] * cyy_n1 * ( dP0dXi[1][ijk_n1] )
                                ) / gbdyn;
                }
            }
            break;
          case CONST_RATE:  //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] * cyy_n1 / gbdyn * (-dP0dXi[1][ijk_n1]);
                }
            }
            break;
         }
     }
  }

  if(nz > 1) {
     double dz0   = g.getDz(0);
     double dz02  = dz0 * dz0;
     double dzn   = g.getDz(nz - 2);
     double gbdz0 = g.getBdz(0);
     double gbdzn = g.getBdz(nz - 1);
     for(j = 0; j < ny; j++ ) {
         j2 = g.toJPerm( j );
         for(i = 0; i < nx; i++ ) {
             i2 = g.toIPerm( i );
             for(k = 1; k < nz - 1; k++ ) {
                 k2 = g.toKPerm( k );
                 k2_p1  = k2 + 1;
                 k2_n1  = k2 - 1;
                 cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                 cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                 ijk_n1 = g.getIndex( i, j, k - 1 );
                 ijk    = g.getIndex( i, j, k     );
                 ijk_p1 = g.getIndex( i, j, k + 1 );
                 permPtr->Perm_z( i, j, k, perm_interface );
                 RHS[ijk] -=   ( perm_interface[2] * cyy_p1 * 
                                  ( dP0dXi[2][ijk   ] + fluidPtr->getGamma(ijk) )
                                -perm_interface[0] * cyy_n1 * 
                                  ( dP0dXi[2][ijk_n1] + fluidPtr->getGamma(ijk) )
                               ) / g.getBdz( k );
             }
         }
     }
     if(g.getGridType() == POINT_DIS) {
        k      = 0 ;
        k2     = g.toKPerm( k );
        k2_p1  = k2 + 1;
        k2_n1  = k2 - 1;
        switch( bType[4] ) {
          case CONST_PRES:
            for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k) ] = 0.0;
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2 = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2   );
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[2]/ dz0 
                                   * ( cyy_p1 * (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                      +bCond[4] * (cyy_p1 - cyy_p0 )/perm_interface[1]
                                     );
                }
            }
            break;
        }
        k      = nz - 1 ;
        k2     = g.toKPerm( k );
        k2_p1  = k2 + 1;
        k2_n1  = k2 - 1;
        switch( bType[5] ) {
          case CONST_PRES:
            for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= 2.* perm_interface[0] / dzn   
                                  * ( cyy_n1 * ( - dP0dXi[2][ijk_n1] - fluidPtr->getGamma(ijk))
                                    - bCond[5]* (cyy_p0 - cyy_n1 )/perm_interface[1]
                                    );
                }
            }
            break;
        }
     }
     else {
        k      = 0 ;
        k2     = g.toKPerm( k );
        k2_p1  = k2 + 1;
        k2_n1  = k2 - 1;
        switch( bType[4] ) {
          case CONST_PRES:  //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * 
                                  (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk) )
                                 -perm_interface[0] * cyy_n1 *
                                  ( (p0[ijk] - bCond[4])/(0.5*gbdz0) + fluidPtr->getGamma(ijk) )
                                ) / gbdz0;
                }
            }
            break;
          case CONST_RATE: //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * 
                                  ( dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                 -perm_interface[0] * cyy_n1 *
                                  ( fluidPtr->getGamma(ijk) - fluidPtr->getGamma(ijk) )
                                ) /gbdz0;
                }
            }
            break;
        }
        k      = nz - 1 ;
        k2     = g.toKPerm( k );
        k2_p1  = k2 + 1;
        k2_n1  = k2 - 1;
        switch( bType[5] ) {
          case CONST_PRES: //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 *
                                  ((bCond[5] - p0[ijk])/(0.5*gbdzn) + fluidPtr->getGamma(ijk))
                                 -perm_interface[0] * cyy_n1 * 
                                  ( dP0dXi[2][ijk_n1] + fluidPtr->getGamma(ijk))
                                 ) / gbdzn;
                }
            }
            break;
          case CONST_RATE: //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[0] * cyy_n1 * 
                                  (-dP0dXi[2][ijk_n1] - fluidPtr->getGamma(ijk))
                                 -perm_interface[2] * cyy_p1 *
                                  ( fluidPtr->getGamma(ijk) - fluidPtr->getGamma(ijk) )
                               ) /gbdzn;
                }
            }
            break;
        }
     }
  }
}     

//============================================================================
void CYPEqn::calcAccu(bool recalA, double dt) {
  /*
  int i, j, k, ijk;
  for(k=0; k<nz; k++) {
      for(j=0; j<ny; j++) {
          for(i=0; i<nx; i++) {
              ijk = i+(j+k*ny)*nx;
              //  before solve, CYP still store the old timestep value
              RHS[ijk] -=  fluidPtr->c[ijk]/dt*CYP[ijkc+ijk*nNode];
              RHS[ijk] -=  permPtr->CY[ijkc+ijk*nNode]
                         *fluidPtr->c[ijk]*dP0dt[ijk];
          }
      }
  }
  */
}

//===== ( output ) =====
void CYPEqn::output( const char* fileName, int flag, const Control &c,
                     int slice_index ) {
  // --- output CYP ----------------------------
  double aa = (unit == FIELD)? psia_pa : 1;
  int i, j, k, ijk;
  int ii = c.getI_Ref();
  int jj = c.getJ_Ref();
  int kk = c.getK_Ref();
  int iijjkk = (ii + (jj+kk*ny)*nx);//[Pipatl]Previous implementation, index based on grid
  //int iijjkk = ((2*ii+1) + (2*jj+1)*nxPerm);//[Pipatl] index based on Y node, only 2D-XY plain
  ofstream os(fileName, ios::out);
  double lcoord;

  os << setiosflags(ios::fixed | ios::showpoint)
     << setprecision(5);

  // (x,y) two dimensional plot
  if( flag == 0 ) {
      for( j = 0; j < ny; j++) {
           for( i = 0; i < nx; i++) {
                ijk = i + j * nx;
                os << g.getX(i) << "  " 
                   << g.getX(j) << "  "
                   << CYP[iijjkk*nNode+ijk] / aa << endl;
           }
          os << endl;
      }
  }

  // (x,z) two dimensional plot
  if( flag == 1 ) {
      for(k = 0; k < nz; k++) {
          for(i = 0; i < nx; i++) {
              ijk = i + k * ny * nx;
              os << g.getX(i) << "  " 
                 << g.getZ(k) << "  "
                 << CYP[iijjkk*nNode+ijk] / aa << endl;
          }
          os << endl;
      }
  }

  // Vertical diagonal surface plot resulting from three dimensional problem.
  // This would work provided nx=ny.
  if( flag == 2 ) {
      for(k = 0; k < nz; k++ ) {
          for(j = 0; j < ny; j++ ) {
              ijk = j+(j+k*ny)*nx;              
              lcoord = sqrt( g.getX(j) * g.getX(j) +
                             g.getY(j) * g.getY(j) );
              os << lcoord << "  " 
                 << g.getZ(k) << "  "
                 << CYP[iijjkk*nNode+ijk] / aa << endl;
          }
          os << endl;
      }
  }

  // Diagonal one dimensional plot resulting from two dimensional problem.
  // This would work provided nx=ny.
  if( flag == 3 ) {
      for(j = 0; j < ny; j++ ) {
          ijk = j + j * nx;
          lcoord = sqrt( g.getX(j) * g.getX(j) +
                         g.getY(j) * g.getY(j) );
          os << lcoord << "  "
             << CYP[iijjkk*nNode+ijk] / aa << endl;
      }
  }

  // Diagonal one dimensional plot resulting from two dimensional problem.
  // This would work provided nx=nz.
  if( flag == 4 ) {
      for(k = 0; k < nz; k++ ) {
          ijk = k + k * nx;
          lcoord = sqrt( g.getX(k) * g.getX(k) +
                         g.getZ(k) * g.getZ(k) );
          os << lcoord << "  "
             << CYP[iijjkk*nNode+ijk] / aa << endl;
      }
  }

  // slice_index belongs to x-direction
  if( flag == 5 ) {
      for(k = 0; k < nz; k++ ) {
          for(j = 0; j < ny; j++ ) {
              ijk = slice_index + (j+k*ny)*nx;
              os << g.getY(j) << "  " 
                 << g.getZ(k) << "  "
                 << CYP[iijjkk*nNode+ijk] / aa << endl;
          }
          os << endl;
      }
  }

  // slice_index belongs to y-direction
  if( flag == 6 ) {
      for(k = 0; k < nz; k++ ) {
          for(i = 0; i < nx; i++ ) {
              ijk = i + (slice_index+k*ny) * nx;
              os << g.getX(i) << "  " 
                 << g.getZ(k) << "  "
                 << CYP[iijjkk * nNode + ijk] / aa << endl;
          }
          os << endl;
      }
  }

  // slice_index belongs to z-direction
  if( flag == 7 ) {
      for(j = 0; j < ny; j++ ) {
          for(i = 0; i < nx; i++ ) {
              ijk = i + (j + slice_index * ny) * nx;
              os << g.getX(i) << "  " 
                 << g.getY(j) << "  "
                 << CYP[iijjkk * nNode + ijk] / aa << endl;
          }
          os << endl;
     }
  }

  // This option is used to save the solution vector in one column.
  if( flag == 8 ) {
      os << nNode << endl;
      for(i = 0; i < nNode; i++ )
          os << CYP[iijjkk * nNode + i] / aa << endl;          
  }
  os.close();
}

//===== ( outputCPY ) =====
//[Pipatl] This function output CYP with pressure at the referrence point and Y everywhere.
void CYPEqn::outputCPY( const char* fileName, int flag, const Control &c,
                     int slice_index ) {
  // --- output CYP ----------------------------
  double aa = (unit == FIELD)? psia_pa : 1;
  int i, j, k, ijk;
  int ii = c.getI_Ref();
  int jj = c.getJ_Ref();
  int kk = c.getK_Ref();
  int iijjkk = (ii + (jj+kk*ny)*nx);//[Pipatl]Previous implementation, index based on grid
  ofstream os(fileName, ios::out);
  double lcoord;

  os << setiosflags(ios::fixed | ios::showpoint)
     << setprecision(5);
  os << "CYP for P at (" << ii+1 << ", " << jj+1 << ", " << kk+1 << ") and Y everywhere." << endl << endl << endl;

  // (x,y) two dimensional plot
  if( flag == 0 ) {
      for( j = 0; j < nyPerm; j++) {
           for( i = 0; i < nxPerm; i++) {
                ijk = i + j * nxPerm;
                os << i << "  " 
                   << j << "  "
                   << CYP[ijk*nNode+iijjkk] / aa << endl;
           }
          os << endl;
      }
  }

  // (x,z) two dimensional plot
  /*if( flag == 1 ) {
      for(k = 0; k < nz; k++) {
          for(i = 0; i < nx; i++) {
              ijk = i + k * ny * nx;
              os << g.getX(i) << "  " 
                 << g.getZ(k) << "  "
                 << CYP[iijjkk*nNode+ijk] / aa << endl;
          }
          os << endl;
      }
  }

  // Vertical diagonal surface plot resulting from three dimensional problem.
  // This would work provided nx=ny.
  if( flag == 2 ) {
      for(k = 0; k < nz; k++ ) {
          for(j = 0; j < ny; j++ ) {
              ijk = j+(j+k*ny)*nx;              
              lcoord = sqrt( g.getX(j) * g.getX(j) +
                             g.getY(j) * g.getY(j) );
              os << lcoord << "  " 
                 << g.getZ(k) << "  "
                 << CYP[iijjkk*nNode+ijk] / aa << endl;
          }
          os << endl;
      }
  }

  // Diagonal one dimensional plot resulting from two dimensional problem.
  // This would work provided nx=ny.
  if( flag == 3 ) {
      for(j = 0; j < ny; j++ ) {
          ijk = j + j * nx;
          lcoord = sqrt( g.getX(j) * g.getX(j) +
                         g.getY(j) * g.getY(j) );
          os << lcoord << "  "
             << CYP[iijjkk*nNode+ijk] / aa << endl;
      }
  }

  // Diagonal one dimensional plot resulting from two dimensional problem.
  // This would work provided nx=nz.
  if( flag == 4 ) {
      for(k = 0; k < nz; k++ ) {
          ijk = k + k * nx;
          lcoord = sqrt( g.getX(k) * g.getX(k) +
                         g.getZ(k) * g.getZ(k) );
          os << lcoord << "  "
             << CYP[iijjkk*nNode+ijk] / aa << endl;
      }
  }

  // slice_index belongs to x-direction
  if( flag == 5 ) {
      for(k = 0; k < nz; k++ ) {
          for(j = 0; j < ny; j++ ) {
              ijk = slice_index + (j+k*ny)*nx;
              os << g.getY(j) << "  " 
                 << g.getZ(k) << "  "
                 << CYP[iijjkk*nNode+ijk] / aa << endl;
          }
          os << endl;
      }
  }

  // slice_index belongs to y-direction
  if( flag == 6 ) {
      for(k = 0; k < nz; k++ ) {
          for(i = 0; i < nx; i++ ) {
              ijk = i + (slice_index+k*ny) * nx;
              os << g.getX(i) << "  " 
                 << g.getZ(k) << "  "
                 << CYP[iijjkk * nNode + ijk] / aa << endl;
          }
          os << endl;
      }
  }

  // slice_index belongs to z-direction
  if( flag == 7 ) {
      for(j = 0; j < ny; j++ ) {
          for(i = 0; i < nx; i++ ) {
              ijk = i + (j + slice_index * ny) * nx;
              os << g.getX(i) << "  " 
                 << g.getY(j) << "  "
                 << CYP[iijjkk * nNode + ijk] / aa << endl;
          }
          os << endl;
     }
  }

  // This option is used to save the solution vector in one column.
  if( flag == 8 ) {
      os << nNode << endl;
      for(i = 0; i < nNode; i++ )
          os << CYP[iijjkk * nNode + i] / aa << endl;          
  }*/
  os.close();
}


//===== printCYP() =====
void CYPEqn::printCYP(){
   ofstream os("CYP.out", ios::out);
   for(int kc = 0; kc < nzPerm; kc++) {
       for(int jc = 0; jc < nyPerm; jc++) {
           for(int ic = 0; ic < nxPerm; ic++) { 
               int ijkc = ic + ( jc + kc * nyPerm ) * nxPerm ;
               for(int i = 0; i < nNode; i++) {
                   os << CYP[i + nNode * ijkc  ] << endl; 
               }
           }
       }
   }
   os.close();
}

void CYPEqn::checkBalance() {
  double perm_interface[3];
  for(int kc = 0; kc < nzPerm; kc++) {
      for(int jc = 0; jc < nyPerm; jc++) {
          for(int ic = 0; ic < nxPerm; ic++) { 
              int ijkc = g.getPermIndex(ic, jc, kc);
	      cout << "ijkc = " << ijkc << endl; 
	      for(int k = 0; k < nz; ++k) {
                  int k2 = g.toKPerm( k ); 
		  int k2_p1 = k2 + 1;
		  int k2_m1 = k2 - 1;
		  double dz = g.getBdz( k );
                  for(int j = 0; j < ny; ++j) {
		      int j2 = g.toJPerm( j );
		      int j2_p1  = j2 + 1;
                      int j2_m1  = j2 - 1;
		      double dy = g.getBdy( j );
                      for(int i = 0; i < nx; i++) {
			  int i2     = g.toIPerm( i );
			  int i2_p1  = i2 + 1;
                          int i2_m1  = i2 - 1;
			  double dx = g.getBdz( i );
			  
	                  int ijk     = g.getIndex(i    ,j    ,k);
	                  int ijk_ip1 = g.getIndex(i + 1,j    ,k);
	                  int ijk_im1 = g.getIndex(i - 1,j    ,k);
	                  int ijk_jp1 = g.getIndex(i    ,j + 1,k);
	                  int ijk_jm1 = g.getIndex(i    ,j - 1,k);
			  int ijk_kp1 = g.getIndex(i    ,j    ,k + 1);
	                  int ijk_km1 = g.getIndex(i    ,j    ,k - 1);

                          double cyp      = CYP[ijk    + ijkc * nNode ];
			  double cyp_ip1  = CYP[ijk_ip1+ ijkc * nNode ];
			  double cyp_im1  = CYP[ijk_im1+ ijkc * nNode ];
			  double cyp_jp1  = CYP[ijk_jp1+ ijkc * nNode ];
			  double cyp_jm1  = CYP[ijk_jm1+ ijkc * nNode ];
			  double cyp_kp1  = CYP[ijk_kp1+ ijkc * nNode ];
			  double cyp_km1  = CYP[ijk_km1+ ijkc * nNode ];

                          double cyy_ip1 = permPtr->getCYY(ic, jc, kc, i2_p1, j2   , k2);
                          double cyy_im1 = permPtr->getCYY(ic, jc, kc, i2_m1, j2   , k2);
                          double cyy_jp1 = permPtr->getCYY(ic, jc, kc, i2   , j2_p1, k2);
                          double cyy_jm1 = permPtr->getCYY(ic, jc, kc, i2   , j2_m1, k2);
                          double cyy_kp1 = permPtr->getCYY(ic, jc, kc, i2   , j2   , k2_p1);
                          double cyy_km1 = permPtr->getCYY(ic, jc, kc, i2   , j2   , k2_m1);
			  
			  permPtr->Perm_x( i, j, k, perm_interface );
			  double Y_ip1  = perm_interface[2];
			  double Y_im1  = perm_interface[0];
			  permPtr->Perm_y( i, j, k, perm_interface );
			  double Y_jp1  = perm_interface[2];
			  double Y_jm1  = perm_interface[0];
			  permPtr->Perm_z( i, j, k, perm_interface );
			  double Y_kp1  = perm_interface[2];
			  double Y_km1  = perm_interface[0];
			  
			  double vx2    = - Y_ip1 * (   ( cyp_ip1 - cyp ) / dx
					              +   cyy_ip1 * dP0dXi[0][ijk ]
						    );
                          double vx1    = - Y_im1 * (   ( cyp - cyp_im1 ) / dx
					              +   cyy_im1 * dP0dXi[0][ijk_im1 ]
					            ); 
	                  double vy2    = - Y_jp1 * (   ( cyp_jp1 - cyp ) / dy
			                              +   cyy_jp1 * dP0dXi[1][ijk ]
				    	            );
                          double vy1    = - Y_jm1 * (   ( cyp - cyp_jm1 ) / dy
					              +   cyy_jm1 * dP0dXi[1][ijk_jm1 ]
					            ); 
			  double vz2    = - Y_kp1 * (   ( cyp_kp1 - cyp ) / dz
			                              +   cyy_kp1 * (  dP0dXi[2][ijk ]
				    	                             + fluidPtr->getGamma(ijk)
						                    )
				                    );
                          double vz1    = - Y_km1 * (   ( cyp - cyp_km1 ) / dz
					              +   cyy_km1 * ( dP0dXi[2][ijk_km1 ]
						                     +fluidPtr->getGamma(ijk_km1)
						                    )
					            ); 
			  if(i == 0     ) vx1 = 0;
			  if(j == 0     ) vy1 = 0;
			  if(k == 0     ) vz1 = 0;
			  if(i == nx - 1) vx2 = 0;
			  if(j == ny - 1) vy2 = 0;
			  if(k == nz - 1) vz2 = 0;

			  double sum = (vx2 - vx1)*dy*dz 
				     + (vy2 - vy1)*dx*dz 
				     + (vz2 - vz1)*dx*dy;
			  if(fabs(sum) > 1.0e-10) {
			     cout << "ijk = " << ijk << ' ' 
				  << "i   = " << i   << ' '
			          << "j   = " << j   << ' '  
				  << sum << endl;
			  }
                      }
                  }
              }
	  }
      }
  }
  exit(0);
}
