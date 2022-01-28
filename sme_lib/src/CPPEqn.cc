/*
 * File: CPPEqn.cc
 * ----------------------------------
 * Implementation for CPPEqn class
 */
#include "CPPEqn.h"
#include "CYPEqn.h"
#include "P0Eqn.h"
#include "Solver.h"
#include "Control.h"
#include "Well.h"

//===== ( CPPEqn() ) =====
CPPEqn::CPPEqn( const P0Eqn &P0ee, const CYPEqn &CYPee, bool debug_ ) 
  : Eqn( P0ee ), P0e( &P0ee ), debug(debug_) {

  if(debug) cout << "CPPEqn::CPPEqn()" <<endl;
  
  p0  = P0e->getP0();
  CYP = CYPee.getCYP();
  CPP         = new double[ nNode * nNode ];
  p_std_Cond  = new double[ nNode ];
  p_std_UnCd  = new double[ nNode ];
  rho_YP_Cond = new double[ nNode ];
  rho_YP_UnCd = new double[ nNode ];
  
  for(int i = 0; i < nNode * nNode; i++) {
      CPP[i] = 0.0;
  }
  for(int i = 0; i < nNode; i++) {
      p_std_Cond[i] = 0.0;
      p_std_UnCd[i] = 0.0;
      rho_YP_Cond[i]= 0.0;
      rho_YP_UnCd[i]= 0.0;
  }
}

//===== (~CPPEqn()) =====
CPPEqn::~CPPEqn() {
  if(debug) cout << "CPPEqn::~CPPEqn()" <<endl;

  delete[] CPP;
  delete[] p_std_Cond;
  delete[] p_std_UnCd;
  delete[] rho_YP_Cond;
  delete[] rho_YP_UnCd;
}

//===== ( Solve ) =====
void CPPEqn::solve( int nWells, const Well *w, const Control &c, 
                    Solver &s, int iter ) {
  // --- currently, no need to record any old value 
  double dt = c.getDt( iter); 
  bool doLU = false;

  // --- calculate CPP cummulatively from time 0 to current time----
  for(int it = 0; it <= iter; it++) {
      // --- point to P0's derivatives for timestep it -------
      if( nx > 1 ) dP0dXi[0] = P0e->getDP0DXi( 0, it ); 
      if( ny > 1 ) dP0dXi[1] = P0e->getDP0DXi( 1, it ); 
      if( nz > 1 ) dP0dXi[2] = P0e->getDP0DXi( 2, it ); 
      if( fluidPtr->isCompressible() ) dP0dt = P0e->getDP0Dt(it);

      // --- update den, poro using P0 of timestep it ----
      if(iter > 0 && (fluidPtr->isCompressible() || 
		      c.getDt( it ) != dt) 
        ) {
         if( fluidPtr->isCompressible() )
         fluidPtr->update( P0e->getP0(it), false );  
         // --- request P0Eqn to recalculate A ---
         ( (P0Eqn*) P0e )->recalcA( c.getBType(), nWells, w, c.getDt( it) );
         dt = c.getDt( it );
         doLU = true;
      }
                      
      for(int kc = 0; kc < nz; kc++) {
          for(int jc = 0; jc < ny; jc++) {
              for(int ic = 0; ic < nx; ic++) {
                  int ijkc = ic + ( jc + kc * ny ) * nx;
                  calcRHS(ijkc, c.getBType(), c.getBCond());
                  if(fluidPtr->isCompressible() )
                     calcAccu( false, c.getDt( it ) );
                  for(int i = 0; i < nNode; i++) 
                      CPP[ijkc*nNode+i] = RHS[i];
              }
          }
      }
      // --- Solve for CPP ----------------
      // --- calculate RHS for each reference point, combine them --- 
      // --- and solve them together by a multi-rhs solver ---------
      s.solve(doLU, *this, nA, CPP);   // X and B share space in CPP
      
      doLU = false ;
   }

   //for(int i = 0; i < nNode; i++) { 
   //    cout << CPP[i * nNode + i] << endl;
   //}
   //exit(0);
   //checkBalance();

   if(debug) printCPP(); 
}

void CPPEqn::calcCondPStd() {
   int ijkc = 0;
   double Y_std = sqrt( permPtr->getYVar(ijkc) );
   for(int i = 0; i < nNode; ++i) {
       if( CPP[i + i * nNode] < 1.e-10 ) {
            p_std_Cond[i] = 0.;
           rho_YP_Cond[i] = 0.;
       } else {
            p_std_Cond[i] = sqrt( CPP[i + i * nNode] );
	   rho_YP_Cond[i] = CYP[i + nNode * ijkc ]/p_std_Cond[i]/Y_std;
       }
   }
}

void CPPEqn::calcUnCdPStd() {
   int ijkc = 0;
   double Y_std = sqrt( permPtr->getYVar(ijkc) );
   for(int i = 0; i < nNode; ++i) {
       if( CPP[i + i * nNode] < 1.e-10 ) {
            p_std_UnCd[i] = 0.;
           rho_YP_UnCd[i] = 0.;
       } else {
            p_std_UnCd[i] = sqrt( CPP[i + i * nNode] );
 	   rho_YP_UnCd[i] = CYP[i + nNode * ijkc ]/p_std_UnCd[i]/Y_std; 
       }
   }
}

//===== (calcRHS() ) =====
void CPPEqn::calcRHS(int& ijkc, const BType* bType, const double *bCond) {
  int nxy = nx * ny;
  int i,  j,  k;
  int i2, j2, k2;
  int ijk, ijk_p1, ijk_n1, ijk_p0;
  double   cyp_p1, cyp_n1, cyp_p0;
  double perm_interface[3];

  for(i = 0; i < nNode; i++) RHS[i] = 0.0;             // set 0

  // --- X direction (dCYPdx) --------------------------
  if( nx > 1) {
      double dx0  = g.getDx( 0 ) ;
      double dx02 = dx0 * dx0 ;
      double dxn  = g.getDx( nx - 2);
      double gbdx0 = g.getBdx( 0 );
      double gbdxn = g.getBdx( nx - 1 );
      for( k = 0; k < nz; k++) {
           k2 = g.toKPerm( k );
           for( j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for( i = 1; i < nx - 1; i++ ) {
                     i2     = g.toIPerm( i );
                     ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                     cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                     ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                     cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                     ijk_n1 =   g.getIndex( i-1, j, k );
                     ijk    =   g.getIndex( i  , j, k );
                     ijk_p1 =   g.getIndex( i+1, j, k );
                     permPtr->Perm_x( i, j, k, perm_interface );
                     RHS[ijk] -= ( perm_interface[2] * cyp_p1 * dP0dXi[0][ijk   ] 
                                  -perm_interface[0] * cyp_n1 * dP0dXi[0][ijk_n1]
                                 ) / g.getBdx( i );
                }
           }
      }
      if(g.getGridType() == POINT_DIS ) {
         i  = 0 ;
         i2 = g.toIPerm( i );
         switch( bType[0] ) {
           case CONST_PRES: 
             for(k = 0; k < nz; k++) 
                for(j = 0; j < ny; j++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
             break;
           case CONST_RATE:
             for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                    ijk_p0 = g.getPermIndex(i2, j2, k2);
                    cyp_p0 = CYP[ ijk_p0 * nNode + ijkc];  
                    ijk    = g.getIndex( i  , j, k );
                    ijk_p1 = g.getIndex( i+1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= 2.*( cyp_p1 * dP0dXi[0][ijk]
                                   + bCond[0] * (cyp_p1 - cyp_p0)/perm_interface[1]
                                   ) / dx0 * perm_interface[2];
                }
             }
             break;
         }
         i  = nx - 1 ;
         //i2     = g.iPerm[i];
         i2 = g.toIPerm( i );
         switch( bType[1] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++) 
                for(j = 0; j < ny; j++) RHS[ g.getIndex(i, j, k ) ] = 0.0;
             break;
           case CONST_RATE:
             for(k = 0; k < nz; k++) {
                //k2 = g.kPerm[k];   
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++) {
                   //j2     = g.jPerm[j];
                   j2     = g.toJPerm( j );
                   ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                   cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                   ijk_p0 = g.getPermIndex(i2, j2, k2);
                   cyp_p0 = CYP[ ijk_p0 * nNode + ijkc];
                   ijk    = g.getIndex(     i, j, k );
                   ijk_n1 = g.getIndex( i - 1, j, k );
                   permPtr->Perm_x( i, j, k, perm_interface );
                   RHS[ijk] -= 2.* perm_interface[0] / dxn   
                                   * ( cyp_n1 * ( -dP0dXi[0][ijk_n1] )
                                   - bCond[1] * (cyp_p0 - cyp_n1 )/ perm_interface[1]
                                   );
                }
             }
             break;
         }
      } else {
         i      = 0 ;
         i2     = g.toIPerm( i );
         switch( bType[0] ) {
           case CONST_PRES: //ok
             for(k = 0; k < nz; k++) {
                 k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyp_p1 * dP0dXi[0][ijk]
                                 -perm_interface[0] * cyp_n1 *
                                    (p0[ijk] - bCond[0])/(0.5*gbdx0)
                                ) / gbdx0;
                }   
             }      
             break;
           case CONST_RATE: //ok
             for(k = 0; k < nz; k++) {
                 k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * cyp_p1 /gbdx0 * dP0dXi[0][ijk];
                }
             }
             break;
         }
         i      = nx - 1;
         i2     = g.toIPerm( i );
         switch( bType[1] ) {
           case CONST_PRES:  //ok
             for(k = 0; k < nz; k++) {
                 k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    ijk_p1 = g.getPermIndex(i2 + 1, j2, k2);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc ];
                    ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyp_p1 *
                                   (bCond[1] - p0[ijk])/(0.5 * gbdxn)
                                 -perm_interface[0] * cyp_n1 * ( dP0dXi[0][ijk_n1])
                                ) / gbdxn;
                }
             }
             break;
           case CONST_RATE: //ok
             for(k = 0; k < nz; k++) {
                 k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    ijk_n1 = g.getPermIndex(i2 - 1, j2, k2);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk    = g.getIndex( i,     j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] * cyp_n1 / gbdxn * (-dP0dXi[0][ijk_n1]);
                }
             }
             break;
         }
      }     
  }

  // --- Y direction (dCYPdy) --------------------------
  if( ny > 1) {
      double dy0  = g.getDy( 0 );
      double dy02 = dy0 * dy0;
      double dyn  = g.getDy( ny - 2 );
      double gbdy0 = g.getBdy( 0 );
      double gbdyn = g.getBdy( ny - 1 );
      for( k = 0; k < nz; k++) {
           k2 = g.toKPerm( k );
           for( i = 0; i < nx; i++) {
                i2     = g.toIPerm( i );
                for( j = 1; j < ny - 1; j++ ) {
                     j2     = g.toJPerm( j );
                     ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                     cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                     ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                     cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                     ijk_n1 =   g.getIndex( i, j-1, k );
                     ijk    =   g.getIndex( i,   j, k );
                     ijk_p1 =   g.getIndex( i, j+1, k );
                     permPtr->Perm_y( i, j, k, perm_interface );
                     RHS[ijk] -= ( perm_interface[2] * cyp_p1 * dP0dXi[1][ijk   ]
                                  -perm_interface[0] * cyp_n1 * dP0dXi[1][ijk_n1]
                                 ) / g.getBdy( j );
                }
           }
      }
      if(g.getGridType() == POINT_DIS ) {
         j      = 0 ;
         //j2     = g.jPerm[j];
         j2 = g.toJPerm( j );
         switch( bType[2] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++)
                for( i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
             break;
           case CONST_RATE:
             for(k = 0; k < nz; k++) {
                //k2 = g.kPerm[k];     
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    //i2     = g.iPerm[i];                        
                    i2     = g.toIPerm( i );
                    ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                    ijk_p0 = g.getPermIndex(i2, j2, k2);
                    cyp_p0 = CYP[ ijk_p0 * nNode + ijkc];
                    ijk    = g.getIndex( i, j  , k );
                    ijk_p1 = g.getIndex( i, j+1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[2] /dy0 
                                   * ( cyp_p1 * dP0dXi[1][ijk]
                                      +bCond[2]* (cyp_p1 - cyp_p0 )/perm_interface[1]
                                     );
                }
             }
             break;
         }
         j      = ny - 1 ;
         //j2     = g.jPerm[j];
         j2 = g.toJPerm( j );
         switch( bType[3] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++)
                for( i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
             break;
           case CONST_RATE:
             for(k = 0; k < nz; k++) {
                //k2 = g.kPerm[k];     
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    //i2     = g.iPerm[i];
                    i2     = g.toIPerm( i );
                    ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk_p0 = g.getPermIndex(i2, j2, k2);
                    cyp_p0 = CYP[ ijk_p0 * nNode + ijkc];
                    ijk    = g.getIndex( i, j     , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[0] / dyn
                                   * ( cyp_n1 * ( -dP0dXi[1][ijk_n1] )
                                     - bCond[3]* (cyp_p0 - cyp_n1)/perm_interface[1]
                                     );
                }
             }
             break;
         }
      } else {
         j      = 0 ;
         j2 = g.toJPerm( j );
         switch( bType[2] ) {
           case CONST_PRES: //ok
             for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                    ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyp_p1 * dP0dXi[1][ijk]
                                 -perm_interface[0] * cyp_n1 *
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
                    ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * cyp_p1 / gbdy0 * dP0dXi[1][ijk];
                }
             }
             break;
         }
         j      = ny - 1;
         //j2     = g.jPerm[j];
         j2 = g.toJPerm( j );
         switch( bType[3] ) {
           case CONST_PRES: //ok
             for(k = 0; k < nz; k++) {
                //k2 = g.kPerm[k];
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    //i2     = g.iPerm[i];
                    i2     = g.toIPerm( i );
                    ijk_p1 = g.getPermIndex(i2, j2 + 1, k2);
                    cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                    ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyp_p1 *
                                    (bCond[3] - p0[ijk])/(0.5*gbdyn)
                                 -perm_interface[0] * cyp_n1 * ( dP0dXi[1][ijk_n1] )
                                ) / gbdyn;
                }
             }
             break;
           case CONST_RATE:  //ok
             for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk_n1 = g.getPermIndex(i2, j2 - 1, k2);
                    cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] * cyp_n1 / gbdyn * (-dP0dXi[1][ijk_n1]);
                }
             }
             break;
         }
      }
  }

  // --- Z direction (dCYPdz, gravity) --------------------------
    if( nz > 1 ) {
      double dz0  = g.getDz( 0 );
      double dz02 = dz0 * dz0;
      double dzn  = g.getDz( nz - 2 );
      double gbdz0 = g.getBdz( 0 );
      double gbdzn = g.getBdz( nz - 1 );
      for( j = 0; j < ny; j++) {
           j2 = g.toJPerm( j );
           for( i = 0; i < nx; i++) {
                i2     = g.toIPerm( i );
                for( k = 1; k < nz - 1; k++ ) {
                     k2 = g.toKPerm( k );
                     ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                     cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                     ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                     cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                     ijk_n1 =   g.getIndex( i, j, k - 1 );
                     ijk    =   g.getIndex( i, j, k     );
                     ijk_p1 =   g.getIndex( i, j, k + 1 );
                     permPtr->Perm_z( i, j, k, perm_interface );
                     RHS[ijk] -=  ( perm_interface[2] * cyp_p1 * 
                                    ( dP0dXi[2][ijk   ] + fluidPtr->getGamma(ijk) )
                                   -perm_interface[0] * cyp_n1 * 
                                    ( dP0dXi[2][ijk_n1] + fluidPtr->getGamma(ijk) )
                                  ) / g.getBdz( k );
                }
           }
      }
      if(g.getGridType() == POINT_DIS) {
         k      = 0 ;
         //k2     = g.kPerm[k];
         k2 = g.toKPerm( k );
         switch( bType[4] ) {
           case CONST_PRES:
             for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k) ] = 0.0;
             break;
           case CONST_RATE:
             for(j = 0; j < ny; j++) {
                //j2 = g.jPerm[j]; 
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                   //i2     = g.iPerm[i]; 
                   i2     = g.toIPerm( i );
                   ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                   cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                   ijk_p0 = g.getPermIndex(i2, j2, k2);
                   cyp_p0 = CYP[ ijk_p0 * nNode + ijkc];
                   ijk    = g.getIndex( i,j,k   );
                   ijk_p1 = g.getIndex( i,j,k+1 );
                   permPtr->Perm_z( i, j, k, perm_interface );
                   RHS[ijk] -= 2. * perm_interface[2] /dz0 
                                  * (  cyp_p1 * (dP0dXi[2][ijk] 
				     + fluidPtr->getGamma( ijk ))
                                    + bCond[4] * (cyp_p1 - cyp_p0 )/perm_interface[1]
                                    );
                }
             }
             break;
         }
         k      = nx - 1 ;
         //k2     = g.kPerm[k];
         k2 = g.toKPerm( k );
         switch( bType[5] ) {
           case CONST_PRES:
             for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k) ] = 0.0;
             break;
           case CONST_RATE:
             for(j = 0; j < ny; j++) {
                //j2 = g.jPerm[j]; 
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                   //i2     = g.iPerm[i]; 
                   i2     = g.toIPerm( i );
                   ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                   cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                   ijk_p0 = g.getPermIndex(i2, j2, k2);
                   cyp_p0 = CYP[ ijk_p0 * nNode + ijkc];
                   ijk    = g.getIndex( i, j, k );
                   ijk_n1 = g.getIndex( i, j, k - 1);
                   permPtr->Perm_z( i, j, k, perm_interface );
                   RHS[ijk] -= 2.* perm_interface[0]/dzn 
                                 *( cyp_n1 * ( - dP0dXi[2][ijk_n1] 
					       - fluidPtr->getGamma(ijk))
                                  - bCond[5]*(cyp_p0 - cyp_n1 )/perm_interface[1]
                                  );
                }
             }
             break;
         }
      } else {
         k      = 0 ;
         //k2     = g.kPerm[k];
         k2 = g.toKPerm( k );
         switch( bType[4] ) {
           case CONST_PRES:  //ok
             for(j = 0; j < ny; j++) {
                //j2 = g.jPerm[j];
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                   //i2     = g.iPerm[i]; 
                   i2     = g.toIPerm( i );
                   ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                   cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                   ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                   cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                   ijk    = g.getIndex( i,j,k   );
                   permPtr->Perm_z( i, j, k, perm_interface );
                   RHS[ijk] -= ( perm_interface[2] * cyp_p1 * 
                                 (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk) )
                                -perm_interface[0] * cyp_n1 *
                                 ( (p0[ijk] - bCond[4])/(0.5*gbdz0) 
				   + fluidPtr->getGamma(ijk) )
                               ) / gbdz0;
                }
             }
             break;
           case CONST_RATE: //ok
             for(j = 0; j < ny; j++) {
                //j2 = g.jPerm[j];
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                  // i2     = g.iPerm[i];
                   i2     = g.toIPerm( i );
                   ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                   cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                   ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                   cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                   ijk    = g.getIndex( i,j,k   );
                   permPtr->Perm_z( i, j, k, perm_interface );
                   RHS[ijk] -= ( perm_interface[2] * cyp_p1 * 
                                 ( dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                -perm_interface[0] * cyp_n1 *
                                 ( fluidPtr->getGamma(ijk) - fluidPtr->getGamma(ijk) )
                               ) /gbdz0;
                }
             }
             break;
         }
         k      = nz - 1 ;
         k2 = g.toKPerm( k );
         switch( bType[5] ) {
           case CONST_PRES: //ok
             for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                   i2     = g.toIPerm( i );
                   ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                   cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                   ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                   cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                   ijk    = g.getIndex( i, j, k );
                   ijk_n1 = g.getIndex( i, j, k - 1);
                   permPtr->Perm_z( i, j, k, perm_interface );
                   RHS[ijk] -= ( perm_interface[2] * cyp_p1 *
                                 (  (bCond[5] - p0[ijk])/(0.5*gbdzn) 
				  + fluidPtr->getGamma(ijk) )
                                -perm_interface[0] * cyp_n1 * 
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
                   ijk_p1 = g.getPermIndex(i2, j2, k2 + 1);
                   cyp_p1 = CYP[ ijk_p1 * nNode + ijkc];
                   ijk_n1 = g.getPermIndex(i2, j2, k2 - 1);
                   cyp_n1 = CYP[ ijk_n1 * nNode + ijkc];
                   ijk    = g.getIndex( i, j, k );
                   ijk_n1 = g.getIndex( i, j, k - 1);
                   permPtr->Perm_z( i, j, k, perm_interface );
                   RHS[ijk] -= ( perm_interface[0] * cyp_n1 * 
                                 (-dP0dXi[2][ijk_n1] - fluidPtr->getGamma(ijk))
                                -perm_interface[2] * cyp_p1 *
                                 ( fluidPtr->getGamma(ijk) - fluidPtr->getGamma(ijk) )
                               ) /gbdzn;
                }
             }
             break;
         }
      }
  }
}

//===== ( calcAccu() ) =====
void CPPEqn::calcAccu(bool recalA, double dt) {
  cout<<"LIYONG: not available!"<<endl;
  /*
  int i, j, k, ijk;
  for(k=0; k<nz; k++) {
  for(j=0; j<ny; j++) {
  for(i=0; i<nx; i++) {
    ijk = i+(j+k*ny)*nx;
    //  before solve, CYP still store the old timestep value
    RHS[ijk] -=  fluidPtr->c[ijk]/dt*CPP[ijk+ijkc*nNode];
    RHS[ijk] -=  CYP[ijk*nNode+ijkc]* fluidPtr->c[ijk]*dP0dt[ijk];
  }}}
  */
}

//=====( output )=====
void CPPEqn::output( const char* fileName, int flag, const Control &c,
                    int slice_index )
{
  int i, j, k, ijk;
  int ii = c.getI_Ref();
  int jj = c.getJ_Ref();
  int kk = c.getK_Ref();
  int iijjkk = (ii+(jj+kk*ny)*nx);
  double aa = CPP[iijjkk+iijjkk*nNode];
  ofstream os(fileName, ios::out);
  double lcoord;

  os << setiosflags(ios::fixed | ios::showpoint)
     << setprecision(5);

  // (x,y) two dimensional plot
  if(flag ==0 ) {
      for(j = 0; j < ny; j++) {
          for(i = 0; i < nx; i++) {
              ijk = i + j * nx;
              if( fabs( aa * CPP[ijk+ijk*nNode]) < 0.0000001 ) 
                os << g.getX( i ) << "  " 
                   << g.getY( j ) << "  " 
		   << 0.0         << "  "
		   << CPP[iijjkk*nNode+ijk]
                   << endl;
              else
                os << g.getX( i ) << "  " 
                   << g.getY( j ) << "  " 
                   << CPP[iijjkk*nNode+ijk] / sqrt( aa*CPP[ijk+ijk*nNode] ) <<"  "
		   << CPP[iijjkk*nNode+ijk]
                   << endl;
          }
          os << endl;
      }
  }

  // (x,z) two dimensional plot
  if(flag == 1 ) {
      for(k = 0; k < nz; k++) {
          for(i = 0; i < nx; i++) {
              ijk = i+k*ny*nx;
              if(aa*CPP[ijk+ijk*nNode]==0 )
                os << g.getX( i ) << "  " 
                   << g.getZ( k ) << "  " 
                   << 0.0 << endl;
              else 
                os << g.getX( i ) << "  " 
                   << g.getZ( k ) << "  "
                   << CPP[iijjkk*nNode+ijk] / sqrt( aa*CPP[ijk+ijk*nNode] )
                   << endl;
          }
          os << endl;
      }
  }

  // Vertical diagonal surface plot resulting from three dimensional problem.
  // This would work provided nx=ny.
  if(flag == 2 ) {
      for( k = 0; k < nz; k++ ) {
          for(j = 0; j < ny; j++ ) {
              ijk = j+(j+k*ny)*nx;              
              lcoord = sqrt(  g.getX( j ) * g.getX( j ) 
                            + g.getY( j ) * g.getY( j ) );
              if( aa*CPP[ijk+ijk*nNode]==0 )
                os << lcoord << "  " 
                   << g.getZ( k ) << "  " 
                   << 0.0 << endl;
              else 
                os << lcoord << "  " 
                   << g.getZ( k ) << "  "
                   << CPP[iijjkk*nNode+ijk] / sqrt( aa*CPP[ijk+ijk*nNode] )
                   << endl;
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
          if( fabs( aa * CPP[ijk+ijk*nNode]) < 0.0000001 )		  
            //os << lcoord << "  " 
	    os << j      << "  "
	       << 0.0    << "  "
	       << CPP[iijjkk*nNode+ijk]
	       << endl;
          else 
            //os << lcoord << "  "
            os << j      << "  "
               << CPP[iijjkk*nNode+ijk] / sqrt( aa*CPP[ijk+ijk*nNode] ) << "  "
	       << CPP[iijjkk*nNode+ijk]
               << endl;
      }
  }

  // Diagonal one dimensional plot resulting from two dimensional problem.
  // This would work provided nx=nz.
  if(flag == 4 ) {
      for(k = 0; k < nz; k++ ) {
          ijk = k+k*nx;
          lcoord = sqrt(  g.getX( k ) * g.getX( k ) 
                        + g.getZ( k ) * g.getZ( k ) );
          if( aa*CPP[ijk+ijk*nNode]==0 )
            os << lcoord << "  " 
               << 0.0 << endl;
          else 
            os << lcoord << "  "
               << CPP[iijjkk*nNode+ijk] / sqrt( aa*CPP[ijk+ijk*nNode] )
               << endl;
      }
  }

  // slice_index belongs to x-direction
  if(flag == 5 ) {
      for(k = 0; k < nz; k++ ) {
          for( j = 0; j < ny; j++ ) {
              ijk = slice_index + (j+k*ny)*nx;
              if( aa*CPP[ijk+ijk*nNode]==0 )
                os << g.getY( j ) << "  " 
                   << g.getZ( k ) << "  " 
                   << 0.0 << endl;
              else 
                os << g.getY( j ) << "  " 
                   << g.getZ( k ) << "  "
                   << CPP[iijjkk*nNode+ijk] / sqrt( aa*CPP[ijk+ijk*nNode] )
                   << endl;
          }
          os << endl;
      }
  }

  // slice_index belongs to y-direction
  if(flag == 6 ) {
      for(k = 0; k < nz; k++ ) {
          for(i = 0; i < nx; i++ ) {
              ijk = i + (slice_index+k*ny)*nx;
              if( aa*CPP[ijk+ijk*nNode]==0 )
                os << g.getX( i ) << "  " 
                   << g.getZ( k ) << "  " 
                   << 0.0 << endl;
              else
                os << g.getX( i ) << "  " 
                   << g.getZ( k ) << "  "
                   << CPP[iijjkk*nNode+ijk] / sqrt( aa*CPP[ijk+ijk*nNode] )
                   << endl;
          }
          os << endl;
      }
  }

  // slice_index belongs to z-direction
  if(flag == 7 ) {
      for(j = 0; j < ny; j++ ) {
          for(i = 0; i < nx; i++ ) {
              ijk = i + (j+slice_index*ny)*nx;
              if( aa*CPP[ijk+ijk*nNode]==0 )
                os << g.getX( i ) << "  " 
                   << g.getY( j ) << "  " 
                   << 0.0 << endl;
              else {
                os << g.getX( i ) << "  " 
                   << g.getY( j ) << "  " 
                   << CPP[iijjkk*nNode+ijk] / sqrt( aa*CPP[ijk+ijk*nNode] )
                   << endl;
              }
          }
          os << endl;
      }
  }

  // This option is used to save the solution vector in one column.
  if(flag == 8 ) {
      os << nNode << endl;
      for(i = 0; i < nNode; i++ ) {
          if( aa*CPP[i+i*nNode]==0 )
            os << 0.0 << endl;
          else
            os << CPP[iijjkk*nNode+i] / sqrt( aa*CPP[i+i*nNode] )
               << endl;
      }          
  }

  os.close();
}

//============================================================================
void CPPEqn::output_PVari( const char* fileName, int flag, int slice_index )
{
  double aa = (unit==FIELD)? psia_pa : 1;
  int i, j, k, ijk;

  ofstream os(fileName, ios::out);
  double lcoord;

  os << setiosflags(ios::fixed | ios::showpoint)
     << setprecision(5);

  // (x,y) two dimensional plot
  if(flag == 0) {
      for(j = 0; j < ny; j++) {
          for(i = 0; i < nx; i++) {
              ijk = i+j*nx;
              os << g.getX( i ) << "  " 
                 << g.getY( j ) << "  "
                 << CPP[ijk*nNode+ijk] / aa / aa << endl;
          }
          os << endl;
      }
  }

  // (x,z) two dimensional plot
  if(flag == 1 ) {
      for(k = 0; k < nz; k++) {
          for(i = 0; i < nx; i++) {
              ijk = i+k*ny*nx;
              os << g.getX( i ) << "  " 
                 << g.getZ( k ) << "  "
                 << CPP[ijk*nNode+ijk] / aa / aa << endl;
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
              os << lcoord << "  " 
                 << g.getZ( k ) << "  "
                 << CPP[ijk*nNode+ijk] / aa / aa << endl;
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
                        + g.getY( j ) * g.getY( j ));
          os << lcoord << "  "
             << CPP[ijk*nNode+ijk] / aa / aa << endl;
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
             << CPP[ijk*nNode+ijk] / aa / aa << endl;
      }
  }

  // slice_index belongs to x-direction
  if(flag == 5 ) {
      for(k = 0; k < nz; k++ ) {
          for( j = 0; j < ny; j++ ) {
              ijk = slice_index + (j+k*ny)*nx;
              os << g.getY( j ) << "  " 
                 << g.getZ( k ) << "  "
                 << CPP[ijk*nNode+ijk] / aa / aa << endl;
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
                 << CPP[ijk*nNode+ijk] / aa / aa << endl;
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
                 << CPP[ijk*nNode+ijk] / aa / aa << endl;
          }
          os << endl;
      }
  }

  // This option is used to save the solution vector in one column.
  if(flag == 8 ) {
      os << nNode << endl;
      for( i=0; i<nNode; i++ )
        os << CPP[i*nNode+i] / aa / aa << endl;          
  }

  os.close();
}

void CPPEqn::printRho(){
   int ijk, ijk1;
   ofstream os("rho.out", ios::out);
   int k1 = nz/2;
   int j1 = ny/2;
   int i1 = nx/2;
   ijk1 = g.getIndex(i1, j1, k1);
   double std1 = sqrt( CPP[ijk1 + ijk1 * nNode]);
   double std;
   for(int k = 0; k < nz; k++) {
       for(int j = 0; j < ny; j++) {
          for(int i = 0; i < nx; i++) {
              ijk = g.getIndex(i, j, k);
              std = sqrt( CPP[ijk + ijk * nNode]);
              os <<i<<' '
                 <<j<<' '
                 <<CPP[ijk + ijk1 * nNode]/std/std1 << endl;
          }
          os<<endl;
       }
   }
   os.close();
}

void CPPEqn::printCPP(){
   ofstream os("CPP.out", ios::out);
   for(int j = 0; j < nNode; j++) 
       os << CPP[j + j * nNode] << endl;
   os<<endl;
   
   for(int j = 0; j < nNode; j++) {
       for(int i = 0; i < nNode; i++) {
           os << CPP[i + j * nNode] << endl;
       }
   }
   os.close();
}

void CPPEqn::checkBalance() {
   double perm_interface[3];
   for(int kc = 0; kc < nz; kc++) {
       for(int jc = 0; jc < ny; jc++) {
          for(int ic = 0; ic < nx; ic++) { 
              int ijkc = g.getIndex(ic, jc, kc);
	      cout << "ijkc = " << ijkc << endl; 
	      for(int k = 0; k < nz; ++k) {
                  int k2    = g.toKPerm( k );
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
			  
			  int ijk2_ip1 = g.getPermIndex(i2_p1,j2,k2);
			  int ijk2_im1 = g.getPermIndex(i2_m1,j2,k2);
			  int ijk2_jp1 = g.getPermIndex(i2,j2_p1,k2);
			  int ijk2_jm1 = g.getPermIndex(i2,j2_m1,k2);
			  int ijk2_kp1 = g.getPermIndex(i2,j2,k2_p1);
			  int ijk2_km1 = g.getPermIndex(i2,j2,k2_m1);

	                  int ijk     = g.getIndex(i    ,j    , k);
	                  int ijk_ip1 = g.getIndex(i + 1,j    , k);
	                  int ijk_im1 = g.getIndex(i - 1,j    , k);
	                  int ijk_jp1 = g.getIndex(i    ,j + 1, k);
	                  int ijk_jm1 = g.getIndex(i    ,j - 1, k);
			  int ijk_kp1 = g.getIndex(i    ,j, k + 1);
	                  int ijk_km1 = g.getIndex(i    ,j, k - 1);
			  
                          double cpp      = CPP[ijk    + ijkc * nNode ];
			  double cpp_ip1  = CPP[ijk_ip1+ ijkc * nNode ];
			  double cpp_im1  = CPP[ijk_im1+ ijkc * nNode ];
			  double cpp_jp1  = CPP[ijk_jp1+ ijkc * nNode ];
			  double cpp_jm1  = CPP[ijk_jm1+ ijkc * nNode ];
			  double cpp_kp1  = CPP[ijk_kp1+ ijkc * nNode ];
			  double cpp_km1  = CPP[ijk_km1+ ijkc * nNode ];
			  
                          double cpy_ip1  = CYP[ijkc   + ijk2_ip1 * nNode ];
                          double cpy_im1  = CYP[ijkc   + ijk2_im1 * nNode ];
                          double cpy_jp1  = CYP[ijkc   + ijk2_jp1 * nNode ];
                          double cpy_jm1  = CYP[ijkc   + ijk2_jm1 * nNode ];
			  double cpy_kp1  = CYP[ijkc   + ijk2_kp1 * nNode ];
                          double cpy_km1  = CYP[ijkc   + ijk2_km1 * nNode ];

			  permPtr->Perm_x( i, j, k, perm_interface );
			  double Y_ip1  = perm_interface[2];
			  double Y_im1  = perm_interface[0];
			  permPtr->Perm_y( i, j, k, perm_interface );
			  double Y_jp1  = perm_interface[2];
			  double Y_jm1  = perm_interface[0];
			  permPtr->Perm_z( i, j, k, perm_interface );
			  double Y_kp1  = perm_interface[2];
			  double Y_km1  = perm_interface[0];

			  double vx2    = - Y_ip1 * (   ( cpp_ip1 - cpp ) / dx
					              +   cpy_ip1 * dP0dXi[0][ijk ]
						    );
                          double vx1    = - Y_im1 * (   ( cpp - cpp_im1 ) / dx
					              +   cpy_im1 * dP0dXi[0][ijk_im1 ]
					            ); 
	                  double vy2    = - Y_jp1 * (   ( cpp_jp1 - cpp ) / dy
			                              +   cpy_jp1 * dP0dXi[1][ijk ]
				    	            );
                          double vy1    = - Y_jm1 * (   ( cpp - cpp_jm1 ) / dy
					              +   cpy_jm1 * dP0dXi[1][ijk_jm1 ]
					            ); 
			  double vz2    = - Y_kp1 * (   ( cpp_kp1 - cpp ) / dz
			                              +   cpy_kp1 * ( dP0dXi[2][ijk ]
							             +fluidPtr->getGamma(ijk_km1)
							            )
				    	            );
                          double vz1    = - Y_km1 * (   ( cpp - cpp_km1 ) / dz
					              +   cpy_km1 * ( dP0dXi[2][ijk_km1 ]
							             +fluidPtr->getGamma(ijk_km1)
							            )
					            ); 

			  if(i == 0) vx1 = 0;
			  if(j == 0) vy1 = 0;
			  if(k == 0) vz1 = 0;
                          if(i == nx - 1) vx2 = 0;
			  if(j == ny - 1) vy2 = 0;
                          if(k == nz - 1) vz2 = 0;
			  double sum = (vx2 - vx1)*dy*dz 
				     + (vy2 - vy1)*dx*dz 
				     + (vz2 - vz1)*dx*dy;

			  //cout << "ijk = " << ijk << ' ' << sum << endl;
			  if(fabs(sum) > 1.0e-15) {
			     cout << "ijk = " << ijk << ' ' 
				  << "i   = " << i   << ' '
			          << "j   = " << j   << ' '  
				  << "k   = " << k   << ' '
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

