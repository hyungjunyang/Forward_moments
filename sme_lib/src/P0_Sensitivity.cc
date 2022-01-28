/*
 * File: P0_Sensitivity.cc
 * ----------------------------------
 * Implementation for P0_Sensitivity class
 */
#include "P0_Sensitivity.h"

//===== P0_Sensitivity( const P0Eqn &P0ee ) =====
P0_Sensitivity::P0_Sensitivity(const P0Eqn &P0ee, MasterPoint *mpts, 
                               PressMeasure *prsm, KrigWeigh *krig)
: Eqn( P0ee ), P0e( &P0ee ), 
  mptsPtr(mpts), prsmPtr(prsm), krigPtr(krig) { 

  debug = false;	  
  if(debug)
     cout << endl
          << "P0_Sensitivity::P0_Sensitivity( )"
	  << endl;
   
  num_meas = prsmPtr->getLength();

  p0  = P0e->getP0();  
  p0_Sensitivity = new double[ num_meas * mptsPtr->getLength() ];
  for(int i = 0; i < num_meas * mptsPtr->getLength(); ++i) {
      p0_Sensitivity[i] = 0.0;
  }
}

//===== ~P0_Sensitivity() =====
P0_Sensitivity::~P0_Sensitivity() {
  if(debug) 
     cout << "P0_Sensitivity::~P0_Sensitivity( )"<<endl;
  delete[] p0_Sensitivity;
}

//===== solve() =====
void P0_Sensitivity::solve(const Control &c, Solver &s) {
  // --- point to P0's derivatives for current timestep -------
  if( nx > 1 ) dP0dXi[0] = P0e->getDP0DXi( 0, 0) ;
  if( ny > 1 ) dP0dXi[1] = P0e->getDP0DXi( 1, 0) ;
  if( nz > 1 ) dP0dXi[2] = P0e->getDP0DXi( 2, 0) ;
  
  // cout << "done at solve () " <<endl;
  for(int m = 0; m < mptsPtr->getLength(); ++m) {
      int ic = mptsPtr->getIm( m );
      int jc = mptsPtr->getJm( m );
      int kc = mptsPtr->getKm( m );
      //cout <<"P0_sen"<< ic <<' '<< jc <<' '<< kc <<' '<<m <<endl;
	     
      calcRHS(m, c.getBType(), c.getBCond());
      s.solve( false, *this, 1, RHS );
      for(int i = 0; i < num_meas; ++i) {
          int ijk = g.getIndex( prsmPtr->getIpmeas(i),
			        prsmPtr->getJpmeas(i),
				prsmPtr->getKpmeas(i)
			      );
          p0_Sensitivity[ i + m * num_meas] = RHS[ ijk ];
      }
  } 
  if(debug) { output();
     cout << "done at solve () " <<endl;
  }
}

//===== ( calcRHS() ) =====
void P0_Sensitivity::calcRHS(int& master, const BType *bType, const double *bCond)
{ 
  int nxy = nx * ny;
  int i , j , k ;
  int i2, j2, k2;
  int ijk, ijk_p1, ijk_n1;
  int i2_p1, i2_n1;
  int j2_p1, j2_n1;
  int k2_p1, k2_n1;
  int      ijk2_p1,   ijk2_n1,   ijk2_p0;
  double   lambda_p1, lambda_n1, lambda_p0;
  double perm_interface[3];
  
  for( i = 0; i < nNode; i++ ) RHS[i] = 0.0;
  
  // --- X direction -------------------------
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
                 i2      = g.toIPerm( i );
                 
                 i2_p1   = i2 + 1;
                 i2_n1   = i2 - 1;
                 ijk2_p1 = g.getPermIndex(i2_p1, j2, k2);
                 ijk2_n1 = g.getPermIndex(i2_n1, j2, k2);
                 lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                 lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
                 
                 ijk     = g.getIndex( i  , j, k );
                 ijk_n1  = g.getIndex( i-1, j, k );
                 permPtr->Perm_x( i, j, k, perm_interface );

                 RHS[ijk] -= ( perm_interface[2] * lambda_p1 * dP0dXi[0][ijk   ]
                              -perm_interface[0] * lambda_n1 * dP0dXi[0][ijk_n1]
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
                    
                    ijk2_p1    = g.getPermIndex(i2_p1, j2, k2);
                    ijk2_p0    = g.getPermIndex(i2   , j2, k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_p0  = krigPtr->getWeigh(master,ijk2_p0);
                    
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= 2.*(  lambda_p1 * dP0dXi[0][ijk]
                                    + bCond[0] * (lambda_p1 - lambda_p0)
                                      /perm_interface[1]
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

                    ijk2_p0    = g.getPermIndex(i2   , j2, k2); 
                      ijk2_n1    = g.getPermIndex(i2_n1, j2, k2);
                    lambda_p0  = krigPtr->getWeigh(master,ijk2_p0);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);

                    ijk    = g.getIndex(     i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= 2.*(  lambda_n1 * ( -dP0dXi[0][ijk_n1] )
                                    - bCond[1] * (lambda_p0 - lambda_n1 )
                                      /perm_interface[1]
                                   ) / dxn * perm_interface[0]; 
                                  
                }
            }
            break;
        }
      }
     else {
        i      = 0 ;
        i2     = g.toIPerm( i );
        i2_p1  = i2 + 1;
        i2_n1  = i2 - 1;
        switch( bType[0] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    
                    ijk2_p1    = g.getPermIndex(i2_p1, j2, k2);
                    ijk2_n1    = g.getPermIndex(i2_n1, j2, k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
                    
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * lambda_p1 * dP0dXi[0][ijk]
                                 -perm_interface[0] * lambda_n1 *
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

                    ijk2_p1    = g.getPermIndex(i2_p1, j2, k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * lambda_p1 / gbdx0 * dP0dXi[0][ijk];
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
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );

                    ijk2_p1    = g.getPermIndex(i2_p1, j2, k2);
                    ijk2_n1    = g.getPermIndex(i2_n1, j2, k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
                    
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * lambda_p1 *
                                   (bCond[1] - p0[ijk])/(0.5 * gbdxn)
                                 -perm_interface[0] * lambda_n1 * ( dP0dXi[0][ijk_n1])
                                ) / gbdxn;
                }
            }
            break;
          case CONST_RATE: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );

                    ijk2_n1    = g.getPermIndex(i2_n1, j2, k2);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
                    
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] * lambda_n1 / gbdxn 
                                * (-dP0dXi[0][ijk_n1]);
                }
            }
            break;
        }
     }
  }

  if(ny > 1 ) {
     double dy0   = g.getDy( 0 );
     double dy02  = dy0 * dy0;
     double dyn   = g.getDy( ny - 2 );
     double gbdy0 = g.getBdy( 0 );
     double gbdyn = g.getBdy( ny - 1 );
     for(k = 0; k < nz; k++ ) {
         k2 = g.toKPerm( k );
         for(j = 1; j < ny - 1; j++ ) {
             j2     = g.toJPerm( j );
             j2_p1  = j2 + 1;
             j2_n1  = j2 - 1;
             for(i = 0; i < nx; i++ ) {
                 i2 = g.toIPerm( i );

                 ijk2_p1    = g.getPermIndex(i2, j2_p1, k2);
                 ijk2_n1    = g.getPermIndex(i2, j2_n1, k2);
                 lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                 lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
                 
                 ijk    = g.getIndex( i, j  , k );
                 ijk_n1 = g.getIndex( i, j-1, k );
                 permPtr->Perm_y( i, j, k, perm_interface );
                 RHS[ijk] -= (  perm_interface[2] * lambda_p1 * dP0dXi[1][ijk    ]
                              - perm_interface[0] * lambda_n1 * dP0dXi[1][ijk_n1 ]
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
                    
                    ijk2_p1    = g.getPermIndex(i2, j2_p1, k2);
                    ijk2_p0    = g.getPermIndex(i2, j2   , k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_p0  = krigPtr->getWeigh(master,ijk2_p0);
                    
                    ijk    = g.getIndex( i, j  , k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[2] /dy0 
                                   * (  lambda_p1 * dP0dXi[1][ijk]
                                      + bCond[2]* (lambda_p1 - lambda_p0 )
                                       /perm_interface[1]
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

                    ijk2_p0    = g.getPermIndex(i2, j2   , k2);
                    ijk2_n1    = g.getPermIndex(i2, j2_n1, k2);
                    lambda_p0  = krigPtr->getWeigh(master,ijk2_p0);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
     
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[0] / dyn
                                   * ( lambda_n1 * ( -dP0dXi[1][ijk_n1] )
                                     - bCond[3]* (lambda_p0 - lambda_n1)
                                       /perm_interface[1]
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

                    ijk2_p1    = g.getPermIndex(i2, j2_p1, k2);
                    ijk2_n1    = g.getPermIndex(i2, j2_n1, k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
    
                    ijk    = g.getIndex( i, j    , k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * lambda_p1 * dP0dXi[1][ijk]
                                 -perm_interface[0] * lambda_n1 *
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
		    
                    ijk2_p1    = g.getPermIndex(i2, j2_p1, k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);

                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * lambda_p1 
			       / gbdy0 * dP0dXi[1][ijk];
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

                    ijk2_p1    = g.getPermIndex(i2, j2_p1, k2);
                    ijk2_n1    = g.getPermIndex(i2, j2_n1, k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
		    
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * lambda_p1 *
                                    (bCond[3] - p0[ijk])/(0.5*gbdyn)
                                 -perm_interface[0] * lambda_n1 * ( dP0dXi[1][ijk_n1] )
                                ) / gbdyn;
                }
            }
            break;
          case CONST_RATE:  //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );

		    ijk2_n1    = g.getPermIndex(i2, j2_n1, k2);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);

                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] * lambda_n1 
			        / gbdyn * (-dP0dXi[1][ijk_n1]);
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
     
     for(k = 1; k < nz - 1; k++ ) {
         k2 = g.toKPerm( k );
         k2_p1  = k2 + 1;
         k2_n1  = k2 - 1;
         for(j = 0; j < ny; j++ ) {
             j2 = g.toJPerm( j );
             for(i = 0; i < nx; i++ ) {
                 i2 = g.toIPerm( i );

                 ijk2_p1    = g.getPermIndex(i2, j2, k2_p1);
                 ijk2_n1    = g.getPermIndex(i2, j2, k2_n1);
                 lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                 lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
		 
                 ijk_n1 = g.getIndex( i, j, k - 1 );
                 ijk    = g.getIndex( i, j, k     );
                 ijk_p1 = g.getIndex( i, j, k + 1 );
                 permPtr->Perm_z( i, j, k, perm_interface );
                 RHS[ijk] -=   ( perm_interface[2] * lambda_p1 * 
                                  ( dP0dXi[2][ijk   ] + fluidPtr->getGamma(ijk) )
                                -perm_interface[0] * lambda_n1 * 
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

                    ijk2_p1    = g.getPermIndex(i2, j2, k2_p1);
                    ijk2_p0    = g.getPermIndex(i2, j2, k2   );
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_p0  = krigPtr->getWeigh(master,ijk2_p0);
		    
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[2]/ dz0 
                                   * ( lambda_p1 * (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                      +bCond[4] * (lambda_p1 - lambda_p0 )/perm_interface[1]
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
		 
                    ijk2_p0    = g.getPermIndex(i2, j2, k2   );
                    ijk2_n1    = g.getPermIndex(i2, j2, k2_n1);
                    lambda_p0  = krigPtr->getWeigh(master,ijk2_p0);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);

                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= 2.* perm_interface[0] / dzn   
                                * ( lambda_n1 * ( - dP0dXi[2][ijk_n1] - fluidPtr->getGamma(ijk))
                                   - bCond[5]* (lambda_p0 - lambda_n1 )/perm_interface[1]
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

                    ijk2_p1    = g.getPermIndex(i2, j2, k2_p1);
                    ijk2_n1    = g.getPermIndex(i2, j2, k2_n1);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
		    
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * lambda_p1 * 
                                  (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk) )
                                 -perm_interface[0] * lambda_n1 *
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

                    ijk2_p1    = g.getPermIndex(i2, j2, k2_p1);
                    ijk2_n1    = g.getPermIndex(i2, j2, k2_n1);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
		    
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * lambda_p1 * 
                                  ( dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                 -perm_interface[0] * lambda_n1 *
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

		    ijk2_p1    = g.getPermIndex(i2, j2, k2_p1);
                    ijk2_n1    = g.getPermIndex(i2, j2, k2_n1);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);

                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);		    
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * lambda_p1 *
                                  ((bCond[5] - p0[ijk])/(0.5*gbdzn) + fluidPtr->getGamma(ijk))
                                 -perm_interface[0] * lambda_n1 * 
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

		    ijk2_p1    = g.getPermIndex(i2, j2, k2_p1);
                    ijk2_n1    = g.getPermIndex(i2, j2, k2_n1);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);

                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[0] * lambda_n1 * 
                                  (-dP0dXi[2][ijk_n1] - fluidPtr->getGamma(ijk))
                                 -perm_interface[2] * lambda_p1 *
                                  ( fluidPtr->getGamma(ijk) - fluidPtr->getGamma(ijk) )
                               ) /gbdzn;
                }
            }
            break;
        }
     }
  }
}     

void P0_Sensitivity::calcAccu(bool recalA, double dt) {
  cout<<"no implementation!" << endl;
}

void P0_Sensitivity::display(){
  for(int m = 0; m < mptsPtr->getLength(); ++m) {
      cout << "m = " << m << ' ';
      for(int i = 0; i < num_meas; ++i) {
          cout << p0_Sensitivity[ i + m * num_meas] <<' ';
      }
      cout <<endl;
  } 
}

void P0_Sensitivity::output() {
  ofstream os("P0_sensitivity.out", ios::out);
//ofstream os("P0_sensitivity.out", ios::app);
  os << setiosflags(ios::fixed | ios::showpoint)
     << setprecision(4);
  for(int m = 0; m < mptsPtr->getLength(); ++m) {
      os << "m = " << m << ' ';
      for(int i = 0; i < num_meas; ++i) {
          os << p0_Sensitivity[ i + m * num_meas] <<' ';
      }
      os <<endl;
  } 
  os.close();
}
