/*
 * File: Velocity.cc
 * ----------------------------------
 * Implementation for Velocity class
 */

#include "Velocity.h"

/* --- Public Methods ----------------------------------- */
// --- Constructors and destructor ---------------
Velocity::Velocity(Grid* grid, Region *regn, Perm *perm, Fluid *fluid, 
                   const P0Eqn &P0ee, const CYPEqn &CYPee, 
                   const CPPEqn &CPPee, const P2Eqn &P2ee, 
                   Unit u, char *Dir, bool debug_) 
 :gridPtr(grid), regnPtr(regn), permPtr(perm), fluidPtr(fluid),
  P0e(&P0ee), debug(debug_) {
  debug = true;
  if(debug) cout << "Velocity::Velocity()" << endl;
  unit  = u; 
  directory = Dir;

  nx    = gridPtr->getNx(); 
  ny    = gridPtr->getNy();
  nz    = gridPtr->getNz();
  nNode = gridPtr->getnNode();

  CYP   = CYPee.getCYP();   
  CPP   = CPPee.getCPP();

  initialize();
}

Velocity::~Velocity() {
  if(debug) cout << "Velocity::~Velocity()" << endl;

  if(nx > 1) {
     delete [] v0[0];
     delete [] Cv1v1;
  }
  if(ny > 1) {
     delete [] v0[1];
     delete [] Cv2v2;
  }
  if(nz > 1) {
     delete [] v0[2];
     delete [] Cv3v3;
  }
  if(nx > 1 && ny > 1) delete [] Cv1v2;
  if(nx > 1 && nz > 1) delete [] Cv1v3; 
  if(ny > 1 && nz > 1) delete [] Cv2v3;
}

//===== initialze data (v0 and Cvivj) =====
void Velocity::initialize() {
  if(nx > 1) {
     v0[0] = new double[ nNode ];
     Cv1v1 = new double[ nNode * nNode ];
  }  else {
     v0[0] = NULL;
     Cv1v1 = NULL;
  }
  if(ny > 1) {
     v0[1] = new double[ nNode ];
     Cv2v2 = new double[ nNode * nNode ];
  } else {
     v0[1] = NULL;
     Cv2v2 = NULL;
  }
  if(nz > 1) {
    v0[2] = new double[ nNode ];
    Cv3v3 = new double[ nNode * nNode ];
  } else {
    v0[2] = NULL;
    Cv3v3 = NULL;    
  }
  if(nx > 1 && ny > 1) {
    Cv1v2 = new double[ nNode * nNode ]; 
  } else {
    Cv1v2 = NULL;
  }
  if(nx > 1 && nz > 1) {
    Cv1v3 = new double[ nNode * nNode ];
  } else {
    Cv1v3 = NULL;
  }
  if(ny > 1 && nz > 1) {
    Cv2v3 = new double[ nNode * nNode ];
  } else {
    Cv2v3 = NULL;
  }
}

void Velocity::solve(const Control &con, int iter) {

  // --- point to P0's derivatives for current timestep ---
  if(nx > 1) dP0dXi[0] = P0e->getDP0DXi(0, iter);
  if(ny > 1) dP0dXi[1] = P0e->getDP0DXi(1, iter);
  if(nz > 1) dP0dXi[2] = P0e->getDP0DXi(2, iter);
   
  // --- calculate v0[] and v2[] ------------------------------
  if(debug) cout << " calculation Mean Velocity" << endl;
  calcVs();
  
  // --- calculate all Cv ------------------------------------
  if(debug) cout << " calculation Velocity Covariances " <<endl;
  calcCViVj(Cv1v1, Cv2v2, Cv1v2);
  
  //check();
  // if(debug) 
  check3d();
  //if(debug) cout << " ok" << endl;
  //
  //if(unit == FIELD) {
  //   cout << "velocity() -> convert()" << endl;
  //   convertUnit();
  //}
}

// ===== liyl ( partially done) ======
// only for cell-center grids and no-flow boundary
void Velocity::calcVs() {
  int nxy = nx * ny;
  double perm_interface[3];
  if(nx > 1) {
     for(int k = 0; k < nz; k++) {
         for(int j = 0; j < ny; j++) {
             for(int i = 0; i < nx - 1; i++) {
                 int ijk = gridPtr->getIndex( i, j, k );
                 permPtr->Perm_x( i, j, k, perm_interface );
                 v0[0][ijk] = - perm_interface[2] * dP0dXi[0][ijk]
                              / fluidPtr->getPoro(ijk);
		                 //(Pipatl) update to second order accuracy (start)(1+\sigma^2/2)
				 int ijkperm = gridPtr->getPermIndex(gridPtr->toIPerm(i)+1,gridPtr->toJPerm(j),gridPtr->toKPerm(k));
				 int ijkp = gridPtr->getIndex( i+1, j, k );
				 v0[0][ijk] -= perm_interface[2]/fluidPtr->getPoro(ijk)*
					 ((0.5*dP0dXi[0][ijk]*permPtr->getCYY(ijkperm,ijkperm))
					 +((CYP[ijkperm*nNode+ijkp]-CYP[ijkperm*nNode+ijk])/gridPtr->getDx(i)));
				 //(Pipatl) update to second order accuracy (end)
			         /*//(Pipatl) update to second order accuracy (start)(exp(\sigma^2/2))
				 int ijkperm = gridPtr->getPermIndex(gridPtr->toIPerm(i)+1,gridPtr->toJPerm(j),gridPtr->toKPerm(k));
				 int ijkp = gridPtr->getIndex( i+1, j, k );
				 v0[0][ijk] -= perm_interface[2]/fluidPtr->getPoro(ijk)*
					 (((exp(0.5*permPtr->getCYY(ijkperm,ijkperm))-1)*dP0dXi[0][ijk])
					 +((CYP[ijkperm*nNode+ijkp]-CYP[ijkperm*nNode+ijk])/gridPtr->getDx(i)));
				 //(Pipatl) update to second order accuracy (end)
		                 //(Pipatl) update to second order accuracy (start)(mixed)
				 int ijkperm = gridPtr->getPermIndex(gridPtr->toIPerm(i)+1,gridPtr->toJPerm(j),gridPtr->toKPerm(k));
				 int ijkp = gridPtr->getIndex( i+1, j, k );
				 v0[0][ijk] -= perm_interface[2]/fluidPtr->getPoro(ijk)*
					 (((0.5*permPtr->getCYY(ijkperm,ijkperm)+exp(0.5*permPtr->getCYY(ijkperm,ijkperm)))*dP0dXi[0][ijk])
					 +(2*(CYP[ijkperm*nNode+ijkp]-CYP[ijkperm*nNode+ijk])/gridPtr->getDx(i)));
				 v0[0][ijk] *= 0.5;
				 //(Pipatl) update to second order accuracy (end)*/
             }
             int i = nx - 1;
             int ijk = gridPtr->getIndex( i, j, k );
             v0[0][ijk] = 0.; 
         }
     }
  }
  if(ny > 1) {
     for(int k = 0; k < nz; k++) {
         for(int j = 0; j < ny - 1; j++) {
             for(int i = 0; i < nx; i++) {
                 int ijk = gridPtr->getIndex( i, j, k );
                 permPtr->Perm_y( i, j, k, perm_interface );
                 v0[1][ijk] = -perm_interface[2] * dP0dXi[1][ijk]
                               / fluidPtr->getPoro(ijk);
		                 //(Pipatl) update to second order accuracy (start)(1+\sigma^2/2)
				 int ijkperm = gridPtr->getPermIndex(gridPtr->toIPerm(i),gridPtr->toJPerm(j)+1,gridPtr->toKPerm(k));
				 int ijkp = gridPtr->getIndex( i, j+1, k );
				 v0[1][ijk] -= perm_interface[2]/fluidPtr->getPoro(ijk)*
					 ((0.5*dP0dXi[1][ijk]*permPtr->getCYY(ijkperm,ijkperm))
					 +((CYP[ijkperm*nNode+ijkp]-CYP[ijkperm*nNode+ijk])/gridPtr->getDy(j)));
				 //(Pipatl) update to second order accuracy (end)
			         /*//(Pipatl) update to second order accuracy (start)(exp(\sigma^2/2))
				 int ijkperm = gridPtr->getPermIndex(gridPtr->toIPerm(i),gridPtr->toJPerm(j)+1,gridPtr->toKPerm(k));
				 int ijkp = gridPtr->getIndex( i, j+1, k );
				 v0[1][ijk] -= perm_interface[2]/fluidPtr->getPoro(ijk)*
				         (((exp(0.5*permPtr->getCYY(ijkperm,ijkperm))-1)*dP0dXi[1][ijk])
					 +((CYP[ijkperm*nNode+ijkp]-CYP[ijkperm*nNode+ijk])/gridPtr->getDy(j)));
				 //(Pipatl) update to second order accuracy (end)
		                 //(Pipatl) update to second order accuracy (start)(mixed)
				 int ijkperm = gridPtr->getPermIndex(gridPtr->toIPerm(i),gridPtr->toJPerm(j)+1,gridPtr->toKPerm(k));
				 int ijkp = gridPtr->getIndex( i, j+1, k );
				 v0[1][ijk] -= perm_interface[2]/fluidPtr->getPoro(ijk)*
					 (((0.5*permPtr->getCYY(ijkperm,ijkperm)+exp(0.5*permPtr->getCYY(ijkperm,ijkperm)))*dP0dXi[1][ijk])
					 +(2*(CYP[ijkperm*nNode+ijkp]-CYP[ijkperm*nNode+ijk])/gridPtr->getDy(j)));
				 v0[1][ijk] *= 0.5;
				 //(Pipatl) update to second order accuracy (end)*/
             }
         }
         int j = ny - 1;
         for(int i = 0; i < nx; i++) {
             int ijk = gridPtr->getIndex( i, j, k );
             v0[1][ijk] = 0.0;
         }
     }
  }
  
  if(nz > 1) {
     for(int k = 0; k < nz - 1; k++) {
         for(int j = 0; j < ny; j++) {
             for(int i = 0; i < nx; i++) {
                 int ijk = gridPtr->getIndex( i, j, k );
                 permPtr->Perm_z( i, j, k, perm_interface );
                 v0[2][ijk] = -(  perm_interface[2] * dP0dXi[2][ijk]
                                + perm_interface[2] * fluidPtr->getGamma(ijk) 
                               ) / fluidPtr->getPoro(ijk);
		                 //(Pipatl) update to second order accuracy (start)(1+\sigma^2/2)
				 int ijkperm = gridPtr->getPermIndex(gridPtr->toIPerm(i),gridPtr->toJPerm(j),gridPtr->toKPerm(k)+1);
				 int ijkp = gridPtr->getIndex( i, j, k+1 );
				 v0[2][ijk] -= perm_interface[2]/fluidPtr->getPoro(ijk)*
					 ((0.5*dP0dXi[2][ijk]*permPtr->getCYY(ijkperm,ijkperm))
					 +((CYP[ijkperm*nNode+ijkp]-CYP[ijkperm*nNode+ijk])/gridPtr->getDz(k)));
				 //(Pipatl) update to second order accuracy (end)
		                 /*//(Pipatl) update to second order accuracy (start)(exp(\sigma^2/2))
				 int ijkperm = gridPtr->getPermIndex(gridPtr->toIPerm(i),gridPtr->toJPerm(j),gridPtr->toKPerm(k)+1);
				 int ijkp = gridPtr->getIndex( i, j, k+1 );
				 v0[2][ijk] -= perm_interface[2]/fluidPtr->getPoro(ijk)*
					 (((exp(0.5*permPtr->getCYY(ijkperm,ijkperm))-1)*dP0dXi[2][ijk])
					 +((CYP[ijkperm*nNode+ijkp]-CYP[ijkperm*nNode+ijk])/gridPtr->getDz(k)));
				 //(Pipatl) update to second order accuracy (end)
		                 //(Pipatl) update to second order accuracy (start)(mixed)
				 int ijkperm = gridPtr->getPermIndex(gridPtr->toIPerm(i),gridPtr->toJPerm(j),gridPtr->toKPerm(k)+1);
				 int ijkp = gridPtr->getIndex( i, j, k+1 );
				 v0[2][ijk] -= perm_interface[2]/fluidPtr->getPoro(ijk)*
					 (((0.5*permPtr->getCYY(ijkperm,ijkperm)+exp(0.5*permPtr->getCYY(ijkperm,ijkperm)))*dP0dXi[2][ijk])
					 +(2*(CYP[ijkperm*nNode+ijkp]-CYP[ijkperm*nNode+ijk])/gridPtr->getDz(k)));
				 v0[2][ijk] *= 0.5;
				 //(Pipatl) update to second order accuracy (end)*/
             }
         }
     }
     int k = nz - 1;
     for(int j = 0; j < ny; j++) {
         for(int i = 0; i < nx; i++) {
             int ijk = gridPtr->getIndex( i, j, k );
             permPtr->Perm_z( i, j, k, perm_interface );
             v0[2][ijk] = 0.0;
         }
     }
  } 
}

// --- calculate liyong li
// notation:
// 1 : point 1
// 2 : point 2
// p : +1

void Velocity::calcCViVj(double *Cv1v1, double *Cv2v2, double *Cv1v2) {
  double perm_1[3],perm_2[3];
  int i1, j1, k1, ijk1, ijk1p;
  int i2, j2, k2, ijk2, ijk2p;
  double v0_1, v0_2;
  double dx, dx1, dx2, dy, dy1, dy2, dz, dz1, dz2;
  double dcpp;
  double dcyp1, dcyp2;
  double cyy;
  int ij0, ij1, ij2, ij3;
  int iPerm1, jPerm1, kPerm1, ijkPerm1, ijkPerm1p;
  int iPerm2, jPerm2, kPerm2, ijkPerm2, ijkPerm2p; 
  int iPerm1p, jPerm1p, kPerm1p;
  int iPerm2p, jPerm2p, kPerm2p;
  double poro1, poro2;
  if(nx > 1) {
     for(k2 = 0; k2 < nz; k2++) {
         kPerm2 = gridPtr->toKPerm(k2);
         for(j2 = 0; j2 < ny; j2++) {
             jPerm2 = gridPtr->toJPerm(j2);
             for(i2 = 0; i2 < nx; i2++) {
                 iPerm2    = gridPtr->toIPerm(i2);
                 iPerm2p   = iPerm2 + 1;
                 ijkPerm2  = gridPtr->getPermIndex(iPerm2    , jPerm2, kPerm2);
                 ijkPerm2p = gridPtr->getPermIndex(iPerm2 + 1, jPerm2, kPerm2);
                 ijk2      = gridPtr->getIndex( i2    , j2, k2 );
                 ijk2p     = gridPtr->getIndex( i2 + 1, j2, k2 );
		 if(i2==nx-1) dx2 = gridPtr->getDx(i2-1);
		 else dx2       = gridPtr->getDx( i2 );
                 v0_2      = v0[0][ijk2];
                 poro2     = fluidPtr->getPoro(ijk2);
                 permPtr->Perm_x( i2, j2, k2, perm_2 );
                 for(k1 = 0; k1 < nz; k1++) {
                     kPerm1 = gridPtr->toKPerm(k1);
                     for(j1 = 0; j1 < ny; j1++) {
                         jPerm1 = gridPtr->toJPerm(j1);
                         for(i1 = 0; i1 < nx; i1++) {
                             iPerm1    = gridPtr->toIPerm(i1);
                             iPerm1p   = iPerm1 + 1;
                             ijkPerm1  = gridPtr->getPermIndex(iPerm1, jPerm1, kPerm1);
                             ijkPerm1p = gridPtr->getPermIndex(iPerm1+1,jPerm1,kPerm1);
                             ijk1      = gridPtr->getIndex( i1    , j1, k1 );
                             ijk1p     = gridPtr->getIndex( i1 + 1, j1, k1 );
			     if(i1==nx-1) dx1 = gridPtr->getDx(i1-1);
			     else dx1  = gridPtr->getDx( i1 );
                             v0_1      = v0[0][ijk1];
                             poro1     = fluidPtr->getPoro(ijk1);
                             permPtr->Perm_x( i1, j1, k1, perm_1);
                             ij0 = ijk1   + ijk2  * nNode;
                             ij3 = ijk1p  + ijk2p * nNode;
                             ij1 = ijk1p  + ijk2  * nNode;
                             ij2 = ijk1   + ijk2p * nNode;
                             if(i1 == nx - 1 || i2 == nx - 1) {
                                Cv1v1[ij0] = 0;
                             } else {
                                dcpp  = (CPP[ij0]+CPP[ij3]-CPP[ij1]-CPP[ij2])/dx1/dx2;
                                dcyp1 = ( CYP[ijkPerm2p * nNode + ijk1p]
                                        - CYP[ijkPerm2p * nNode + ijk1 ])/dx1;
                                dcyp2 = ( CYP[ijkPerm1p * nNode + ijk2p]
                                        - CYP[ijkPerm1p * nNode + ijk2 ])/dx2;
                                cyy = permPtr->getCYY(iPerm1p, jPerm1, kPerm1, 
                                                      iPerm2p, jPerm2, kPerm2);
                                /*Cv1v1[ij0] = dcpp * perm_1[2] * perm_2[2] / poro1/poro2
                                           + cyy * v0_1 * v0_2
                                           - perm_1[2]  * v0_2 * dcyp1 / poro1
                                           - perm_2[2]  * v0_1 * dcyp2 / poro2;*/
								//(Pipatl) update to second order accuracy (start)
								double dcyp11 = ( CYP[ijkPerm1p * nNode + ijk1p]
												-CYP[ijkPerm1p * nNode + ijk1 ])/dx1;
								double dcyp22 = ( CYP[ijkPerm2p * nNode + ijk2p]
												-CYP[ijkPerm2p * nNode + ijk2 ])/dx2;
								Cv1v1[ij0] = perm_1[2]*perm_2[2]/poro1/poro2*
									(
										(cyy*dP0dXi[0][ijk1]*dP0dXi[0][ijk2])+
										(dP0dXi[0][ijk1]*dcyp2)+
										(dP0dXi[0][ijk2]*dcyp1)+
										dcpp/*-
										(	0.25*
											((dP0dXi[0][ijk1]*permPtr->getCYY(ijkPerm1p,ijkPerm1p))+(2*dcyp11))*
											((dP0dXi[0][ijk2]*permPtr->getCYY(ijkPerm2p,ijkPerm2p))+(2*dcyp22))
										)*/
									);
//cout<<(cyy*dP0dXi[0][ijk1]*dP0dXi[0][ijk2])<<" "<<(dP0dXi[0][ijk1]*dcyp2)<<" "<<(dP0dXi[0][ijk2]*dcyp1)<<" "<<dcpp<<" "<<0.25*((dP0dXi[0][ijk1]*permPtr->getCYY(ijkPerm1p,ijkPerm1p))+(2*dcyp11))*((dP0dXi[0][ijk2]*permPtr->getCYY(ijkPerm2p,ijkPerm2p))+(2*dcyp22));
//cout<<endl;
								//(Pipatl) update to second order accuracy (end)*/
                             }
                         }
                      }
                  }
              }
          }
      }
  }

  if(ny > 1) {
     for(k2 = 0; k2 < nz; k2++) {
         kPerm2 = gridPtr->toKPerm(k2);
         for(j2 = 0; j2 < ny; j2++) {
             jPerm2  = gridPtr->toJPerm(j2);
             jPerm2p = jPerm2  + 1;
	     if(j2==ny-1) dy2 = gridPtr->getDy(j2-1);
             else dy2 = gridPtr->getDy( j2 );
             for(i2 = 0; i2 < nx; i2++) {
                 iPerm2    = gridPtr->toIPerm(i2);
                 ijkPerm2  = gridPtr->getPermIndex( iPerm2, jPerm2    , kPerm2);
                 ijkPerm2p = gridPtr->getPermIndex( iPerm2, jPerm2 + 1, kPerm2);
                 ijk2      = gridPtr->getIndex( i2, j2    , k2 );
                 ijk2p     = gridPtr->getIndex( i2, j2 + 1, k2 );
                 v0_2      = v0[1][ijk2];
                 poro2     = fluidPtr->getPoro(ijk2);
                 permPtr->Perm_y( i2, j2, k2, perm_2 );
                 for(k1 = 0; k1 < nz; k1++) {
                     kPerm1 = gridPtr->toKPerm(k1);
                     for(j1 = 0; j1 < ny; j1++) {
                         jPerm1  = gridPtr->toJPerm(j1);
                         jPerm1p = jPerm1 + 1;
			 if(j1==ny-1) dy1 = gridPtr->getDy(j1-1);
                         else dy1 = gridPtr->getDy( j1 );
                         for(i1 = 0; i1 < nx; i1++) {
                             iPerm1    = gridPtr->toIPerm(i1);
                             ijkPerm1  = gridPtr->getPermIndex(iPerm1,jPerm1, kPerm1);
                             ijkPerm1p = gridPtr->getPermIndex(iPerm1,jPerm1+1,kPerm1);
                             ijk1      = gridPtr->getIndex( i1, j1    , k1 );
                             ijk1p     = gridPtr->getIndex( i1, j1 + 1, k1 );
                             v0_1      = v0[1][ijk1];
                             poro1     = fluidPtr->getPoro(ijk1);
                             permPtr->Perm_y( i1, j1, k1, perm_1);
                             ij0 = ijk1   + ijk2  * nNode;
                             ij3 = ijk1p  + ijk2p * nNode;
                             ij1 = ijk1p  + ijk2  * nNode;
                             ij2 = ijk1   + ijk2p * nNode;
                             if(j1 == ny - 1 || j2 == ny - 1) {
                                Cv2v2[ij0] = 0.;     
                             } else {
                                dcpp  = ( CPP[ij0]+CPP[ij3]-CPP[ij1]-CPP[ij2])/dy2/dy1;
                                dcyp1 = ( CYP[ijkPerm2p * nNode + ijk1p]
                                        - CYP[ijkPerm2p * nNode + ijk1 ])/dy1;
                                dcyp2 = ( CYP[ijkPerm1p * nNode + ijk2p]
                                        - CYP[ijkPerm1p * nNode + ijk2 ])/dy2;
                                cyy = permPtr->getCYY(iPerm1, jPerm1p, kPerm1, 
                                                      iPerm2, jPerm2p, kPerm2);
                                /*Cv2v2[ij0] = dcpp * perm_1[2] * perm_2[2]/poro1/poro2 
                                           + cyy * v0_1 * v0_2
                                           - perm_1[2] * v0_2 * dcyp1/poro1
                                           - perm_2[2] * v0_1 * dcyp2/poro2;*/
								//(Pipatl) update to second order accuracy (start)
								double dcyp11 = ( CYP[ijkPerm1p * nNode + ijk1p]
												-CYP[ijkPerm1p * nNode + ijk1 ])/dy1;
								double dcyp22 = ( CYP[ijkPerm2p * nNode + ijk2p]
												-CYP[ijkPerm2p * nNode + ijk2 ])/dy2;
								Cv2v2[ij0] = perm_1[2]*perm_2[2]/poro1/poro2*
									(
										(cyy*dP0dXi[1][ijk1]*dP0dXi[1][ijk2])+
										(dP0dXi[1][ijk1]*dcyp2)+
										(dP0dXi[1][ijk2]*dcyp1)+
										dcpp/*-
										(	0.25*
											((dP0dXi[1][ijk1]*permPtr->getCYY(ijkPerm1p,ijkPerm1p))+(2*dcyp11))*
											((dP0dXi[1][ijk2]*permPtr->getCYY(ijkPerm2p,ijkPerm2p))+(2*dcyp22))
										)*/
									);
								//(Pipatl) update to second order accuracy (end)
                             }
                         }
                      }
                  }
              }
          }
      }
  }

  if(nz > 1) {
     for(k2 = 0; k2 < nz; k2++) {
         kPerm2 = gridPtr->toKPerm(k2);
         kPerm2p= kPerm2 + 1;
         dz2    = gridPtr->getDz( k2 );
         for(j2 = 0; j2 < ny; j2++) {
             jPerm2 = gridPtr->toJPerm(j2);
             for(i2 = 0; i2 < nx; i2++) {
                 iPerm2    = gridPtr->toIPerm(i2);
                 ijkPerm2  = gridPtr->getPermIndex( iPerm2, jPerm2, kPerm2    );
                 ijkPerm2p = gridPtr->getPermIndex( iPerm2, jPerm2, kPerm2 + 1);
                 ijk2      = gridPtr->getIndex( i2, j2, k2    );
                 ijk2p     = gridPtr->getIndex( i2, j2, k2 + 1);
                 v0_2      = v0[2][ijk2];
                 poro2     = fluidPtr->getPoro(ijk2);
                 permPtr->Perm_z( i2, j2, k2, perm_2 );
                 for(k1 = 0; k1 < nz; k1++) {
                     kPerm1 = gridPtr->toKPerm(k1);
                     kPerm1p= kPerm1 + 1;
                     dz1    = gridPtr->getDz( k1 );
                     for(j1 = 0; j1 < ny; j1++) {
                         jPerm1 = gridPtr->toJPerm(j1);
                         for(i1 = 0; i1 < nx; i1++) {
                             iPerm1    = gridPtr->toIPerm(i1);
                             ijkPerm1  = gridPtr->getPermIndex(iPerm1,jPerm1,kPerm1    );
                             ijkPerm1p = gridPtr->getPermIndex(iPerm1,jPerm1,kPerm1 + 1);
                             ijk1      = gridPtr->getIndex( i1, j1, k1    );
                             ijk1p     = gridPtr->getIndex( i1, j1, k1 + 1);
                             v0_1      = v0[2][ijk1];
                             poro1     = fluidPtr->getPoro(ijk1);
                             permPtr->Perm_z( i1, j1, k1, perm_1);
                             ij0 = ijk1   + ijk2  * nNode;
                             ij3 = ijk1p  + ijk2p * nNode;
                             ij1 = ijk1p  + ijk2  * nNode;
                             ij2 = ijk1   + ijk2p * nNode;
                             if(k1 == nz - 1 || k2 == nz - 1) {
                                Cv3v3[ij0] = 0;
                             } else {
                                dcpp  = ( CPP[ij0]+CPP[ij3]-CPP[ij1]-CPP[ij2])/dz2/dz1;
                                dcyp1 = ( CYP[ijkPerm2p * nNode + ijk1p]
                                        - CYP[ijkPerm2p * nNode + ijk1 ])/dz1;
                                dcyp2 = ( CYP[ijkPerm1p * nNode + ijk2p]
                                        - CYP[ijkPerm1p * nNode + ijk2 ])/dz2;
                                cyy = permPtr->getCYY(iPerm1, jPerm1, kPerm1p, 
                                                      iPerm2, jPerm2, kPerm2p );
                                /*Cv3v3[ij0] = dcpp * perm_1[2] * perm_2[2] / poro1 /poro2
                                           + cyy * v0_1 * v0_2
                                           - perm_1[2] * v0_2 * dcyp1/poro1
                                           - perm_2[2] * v0_1 * dcyp2/poro2;*/
								//(Pipatl) update to second order accuracy (start)
								double dcyp11 = ( CYP[ijkPerm1p * nNode + ijk1p]
												-CYP[ijkPerm1p * nNode + ijk1 ])/dz1;
								double dcyp22 = ( CYP[ijkPerm2p * nNode + ijk2p]
												-CYP[ijkPerm2p * nNode + ijk2 ])/dz2;
								Cv3v3[ij0] = perm_1[2]*perm_2[2]/poro1/poro2*
									(
										(cyy*dP0dXi[2][ijk1]*dP0dXi[2][ijk2])+
										(dP0dXi[2][ijk1]*dcyp2)+
										(dP0dXi[2][ijk2]*dcyp1)+
										dcpp/*-
										(	0.25*
											((dP0dXi[2][ijk1]*permPtr->getCYY(ijkPerm1p,ijkPerm1p))+(2*dcyp11))*
											((dP0dXi[2][ijk2]*permPtr->getCYY(ijkPerm2p,ijkPerm2p))+(2*dcyp22))
										)*/
									);
								//(Pipatl) update to second order accuracy (end)
                             }
                         }
                      }
                  }
              }
          }
      }
  }
  //Cv1v2
  if(nx > 1 && ny > 1) {
     for(k2 = 0; k2 < nz; k2++) {
         kPerm2 = gridPtr->toKPerm(k2);
         for(j2 = 0; j2 < ny; j2++) {
             jPerm2  = gridPtr->toJPerm(j2);
             jPerm2p = jPerm2 + 1;
	     if(j2==ny-1) dy2 = gridPtr->getDy(j2-1);
             else dy2 = gridPtr->getDy( j2 );
             for(i2 = 0; i2 < nx; i2++) {
                 iPerm2    = gridPtr->toIPerm(i2);
                 ijkPerm2  = gridPtr->getPermIndex( iPerm2, jPerm2    , kPerm2);
                 ijkPerm2p = gridPtr->getPermIndex( iPerm2, jPerm2 + 1, kPerm2);
                 ijk2      = gridPtr->getIndex( i2, j2    , k2 );
                 ijk2p     = gridPtr->getIndex( i2, j2 + 1, k2 );
                 v0_2      = v0[1][ijk2];
                 poro2     = fluidPtr->getPoro(ijk2);
                 permPtr->Perm_y( i2, j2, k2, perm_2 );
                 for(k1 = 0; k1 < nz; k1++) {
                     kPerm1 = gridPtr->toKPerm(k1);
                     for(j1 = 0; j1 < ny; j1++) {
                         jPerm1 = gridPtr->toJPerm(j1);
                         for(i1 = 0; i1 < nx; i1++) {
                             iPerm1    = gridPtr->toIPerm(i1);
                             iPerm1p   = iPerm1 + 1;
			     if(i1==nx-1) dx1 = gridPtr->getDx(i1-1);
                             else dx1  = gridPtr->getDx( i1 );
                             ijkPerm1  = gridPtr->getPermIndex(iPerm1,jPerm1,kPerm1);
                             ijkPerm1p = gridPtr->getPermIndex(iPerm1+1,jPerm1,kPerm1);
                             ijk1      = gridPtr->getIndex( i1    , j1, k1 );
                             ijk1p     = gridPtr->getIndex( i1 + 1, j1, k1 );
                             v0_1      = v0[0][ijk1];
                             poro1     = fluidPtr->getPoro(ijk1);
                             permPtr->Perm_x( i1, j1, k1, perm_1);
                             ij0 = ijk1   + ijk2  * nNode;
                             ij3 = ijk1p  + ijk2p * nNode;
                             ij1 = ijk1p  + ijk2  * nNode;
                             ij2 = ijk1   + ijk2p * nNode;
                             
                             if(i1 == nx - 1 || j2 == ny - 1) {
                                Cv1v2[ij0] = 0.;
                             } else {
                                dcpp  = (CPP[ij0]+CPP[ij3]-CPP[ij1]-CPP[ij2])/dy2/dx1;
                                dcyp1 = ( CYP[ijkPerm2p * nNode + ijk1p]
                                        - CYP[ijkPerm2p * nNode + ijk1 ])/dx1;
                                dcyp2 = ( CYP[ijkPerm1p * nNode + ijk2p]
                                        - CYP[ijkPerm1p * nNode + ijk2 ])/dy2;
                                cyy = permPtr->getCYY(iPerm1p, jPerm1 , kPerm1, 
                                                      iPerm2 , jPerm2p, kPerm2);
                                /*Cv1v2[ij0] = dcpp * perm_1[2] * perm_2[2]/poro1/poro2
                                           + cyy * v0_1 * v0_2
                                           - perm_1[2] * v0_2 * dcyp1/poro1
                                           - perm_2[2] * v0_1 * dcyp2/poro2;*/
								//(Pipatl) update to second order accuracy (start)
								double dcyp11 = ( CYP[ijkPerm1p * nNode + ijk1p]
												-CYP[ijkPerm1p * nNode + ijk1 ])/dx1;
								double dcyp22 = ( CYP[ijkPerm2p * nNode + ijk2p]
												-CYP[ijkPerm2p * nNode + ijk2 ])/dy2;
								Cv1v2[ij0] = perm_1[2]*perm_2[2]/poro1/poro2*
									(
										(cyy*dP0dXi[0][ijk1]*dP0dXi[1][ijk2])+
										(dP0dXi[0][ijk1]*dcyp2)+
										(dP0dXi[1][ijk2]*dcyp1)+
										dcpp/*-
										(	0.25*
											((dP0dXi[0][ijk1]*permPtr->getCYY(ijkPerm1p,ijkPerm1p))+(2*dcyp11))*
											((dP0dXi[1][ijk2]*permPtr->getCYY(ijkPerm2p,ijkPerm2p))+(2*dcyp22))
										)*/
									);
								//(Pipatl) update to second order accuracy (end)

                             }
                         }
                      }
                  } 
              }
          }
      }
  }
  //cout << "ok2" << endl;
  //Cv1v3
  if(nx > 1 && nz > 1) {
     for(k2 = 0; k2 < nz; k2++) {
         kPerm2 = gridPtr->toKPerm(k2);
         kPerm2p = kPerm2 + 1;
         dz2    = gridPtr->getDz(  k2);
         for(j2 = 0; j2 < ny; j2++) {
             jPerm2 = gridPtr->toJPerm(j2);
             for(i2 = 0; i2 < nx; i2++) {
                 iPerm2    = gridPtr->toIPerm(i2);
                 ijkPerm2  = gridPtr->getPermIndex( iPerm2, jPerm2, kPerm2    );
                 ijkPerm2p = gridPtr->getPermIndex( iPerm2, jPerm2, kPerm2 + 1);
                 ijk2      = gridPtr->getIndex( i2, j2, k2    );
                 ijk2p     = gridPtr->getIndex( i2, j2, k2 + 1);
                 v0_2      = v0[2][ijk2];
                 poro2     = fluidPtr->getPoro(ijk2);
                 permPtr->Perm_z( i2, j2, k2, perm_2 );

                 for(k1 = 0; k1 < nz; k1++) {
                     kPerm1 = gridPtr->toKPerm(k1);
                     for(j1 = 0; j1 < ny; j1++) {
                         jPerm1 = gridPtr->toJPerm(j1);
                         for(i1 = 0; i1 < nx; i1++) {
                             iPerm1    = gridPtr->toIPerm(i1);
                             iPerm1p   = iPerm1 + 1;
                             dx1       = gridPtr->getDx(  i1);
                             ijkPerm1  = gridPtr->getPermIndex(iPerm1  ,jPerm1,kPerm1);
                             ijkPerm1p = gridPtr->getPermIndex(iPerm1+1,jPerm1,kPerm1);
                             ijk1      = gridPtr->getIndex( i1    , j1, k1 );
                             ijk1p     = gridPtr->getIndex( i1 + 1, j1, k1 );
                             v0_1      = v0[0][ijk1];
                             poro1     = fluidPtr->getPoro(ijk1);
                             permPtr->Perm_x( i1, j1, k1, perm_1);

                             ij0 = ijk1   + ijk2  * nNode;
                             ij3 = ijk1p  + ijk2p * nNode;
                             ij1 = ijk1p  + ijk2  * nNode;
                             ij2 = ijk1   + ijk2p * nNode;
                             if(i1 == nx - 1 || k2 == nz - 1) {
                                Cv1v3[ij0] = 0;
                             } else {
                                dcpp  = (CPP[ij0]+CPP[ij3]-CPP[ij1]-CPP[ij2])/dz2/dx1;
                                dcyp1 = ( CYP[ijkPerm2p * nNode + ijk1p]
                                        - CYP[ijkPerm2p * nNode + ijk1 ] )/dx1;
                                dcyp2 = ( CYP[ijkPerm1p * nNode + ijk2p]
                                        - CYP[ijkPerm1p * nNode + ijk2 ] )/dz2;
                                cyy = permPtr->getCYY(iPerm1p, jPerm1, kPerm1, 
                                                      iPerm2 , jPerm2, kPerm2p);
                                /*Cv1v3[ij0] = dcpp * perm_1[2] * perm_2[2] / poro1/poro2
                                           + cyy * v0_1 * v0_2
                                           - perm_1[2] * v0_2 * dcyp1/poro1
                                           - perm_2[2] * v0_1 * dcyp2/poro2;*/
								//(Pipatl) update to second order accuracy (start)
								double dcyp11 = ( CYP[ijkPerm1p * nNode + ijk1p]
												-CYP[ijkPerm1p * nNode + ijk1 ])/dx1;
								double dcyp22 = ( CYP[ijkPerm2p * nNode + ijk2p]
												-CYP[ijkPerm2p * nNode + ijk2 ])/dz2;
								Cv1v3[ij0] = perm_1[2]*perm_2[2]/poro1/poro2*
									(
										(cyy*dP0dXi[0][ijk1]*dP0dXi[2][ijk2])+
										(dP0dXi[0][ijk1]*dcyp2)+
										(dP0dXi[2][ijk2]*dcyp1)+
										dcpp/*-
										(	0.25*
											((dP0dXi[0][ijk1]*permPtr->getCYY(ijkPerm1p,ijkPerm1p))+(2*dcyp11))*
											((dP0dXi[2][ijk2]*permPtr->getCYY(ijkPerm2p,ijkPerm2p))+(2*dcyp22))
										)*/
									);
								//(Pipatl) update to second order accuracy (end)
                             }
                         }
                      }
                  } 
              }
          }
      }
  }
  //cout << " ok3 " << endl;
  //Cv2v3
  if(ny > 1 && nz > 1) {
     for(k2 = 0; k2 < nz; k2++) {
         kPerm2  = gridPtr->toKPerm(k2);
         kPerm2p = kPerm2 + 1;
         dz2    = gridPtr->getDz( k2 );
         for(j2 = 0; j2 < ny; j2++) {
             jPerm2 = gridPtr->toJPerm(j2);
             for(i2 = 0; i2 < nx; i2++) {
                 iPerm2    = gridPtr->toIPerm(i2);
                 ijkPerm2  = gridPtr->getPermIndex( iPerm2, jPerm2, kPerm2    );
                 ijkPerm2p = gridPtr->getPermIndex( iPerm2, jPerm2, kPerm2 + 1);
                 ijk2      = gridPtr->getIndex( i2, j2, k2    );
                 ijk2p     = gridPtr->getIndex( i2, j2, k2 + 1);
                 v0_2      = v0[2][ijk2];
                 poro2     = fluidPtr->getPoro(ijk2);
                 permPtr->Perm_z( i2, j2, k2, perm_2 );

                 for(k1 = 0; k1 < nz; k1++) {
                     kPerm1 = gridPtr->toKPerm(k1);
                     for(j1 = 0; j1 < ny; j1++) {
                         jPerm1  = gridPtr->toJPerm(j1);
                         jPerm1p = jPerm1 + 1;
                         dy1       = gridPtr->getDy( j1 );
                         for(i1 = 0; i1 < nx; i1++) {
                             iPerm1    = gridPtr->toIPerm(i1);
                             
                             ijkPerm1  = gridPtr->getPermIndex(iPerm1,jPerm1  ,kPerm1);
                             ijkPerm1p = gridPtr->getPermIndex(iPerm1,jPerm1+1,kPerm1);
                             ijk1      = gridPtr->getIndex( i1, j1  , k1 );
                             ijk1p     = gridPtr->getIndex( i1, j1+1, k1 );
                             v0_1      = v0[1][ijk1];
                             poro1     = fluidPtr->getPoro(ijk1);
                             permPtr->Perm_y( i1, j1, k1, perm_1);

                             ij0 = ijk1   + ijk2  * nNode;
                             ij3 = ijk1p  + ijk2p * nNode;
                             ij1 = ijk1p  + ijk2  * nNode;
                             ij2 = ijk1   + ijk2p * nNode;
                             if(j1 == ny - 1 || k2 == nz - 1) {
                                Cv2v3[ij0] = 0;
                             } else {
                                dcpp  = (CPP[ij0]+CPP[ij3]-CPP[ij1]-CPP[ij2])/dz2/dy1;
                                dcyp1 = ( CYP[ijkPerm2p * nNode + ijk1p]
                                        - CYP[ijkPerm2p * nNode + ijk1 ])/dy1;
                                dcyp2 = ( CYP[ijkPerm1p * nNode + ijk2p]
                                        - CYP[ijkPerm1p * nNode + ijk2 ])/dz2;
                                cyy = permPtr->getCYY(iPerm1, jPerm1p, kPerm1, 
                                                      iPerm2, jPerm2 , kPerm2p);
                                /*Cv2v3[ij0] = dcpp * perm_1[2] * perm_2[2] /poro1/poro2
                                           + cyy * v0_1 * v0_2
                                           - perm_1[2] * v0_2 * dcyp1/poro1
                                           - perm_2[2] * v0_1 * dcyp2/poro2;*/
								//(Pipatl) update to second order accuracy (start)
								double dcyp11 = ( CYP[ijkPerm1p * nNode + ijk1p]
												-CYP[ijkPerm1p * nNode + ijk1 ])/dy1;
								double dcyp22 = ( CYP[ijkPerm2p * nNode + ijk2p]
												-CYP[ijkPerm2p * nNode + ijk2 ])/dz2;
								Cv2v3[ij0] = perm_1[2]*perm_2[2]/poro1/poro2*
									(
										(cyy*dP0dXi[1][ijk1]*dP0dXi[2][ijk2])+
										(dP0dXi[1][ijk1]*dcyp2)+
										(dP0dXi[2][ijk2]*dcyp1)+
										dcpp/*-
										(	0.25*
											((dP0dXi[1][ijk1]*permPtr->getCYY(ijkPerm1p,ijkPerm1p))+(2*dcyp11))*
											((dP0dXi[2][ijk2]*permPtr->getCYY(ijkPerm2p,ijkPerm2p))+(2*dcyp22))
										)*/
									);
								//(Pipatl) update to second order accuracy (end)
                             }
                         }
                      }
                  } 
              }
          }
      }
  }
}

void Velocity::convertUnit() {
  double aaa  =     ft_m/day_sec ;
  double aaa2 = pow(ft_m/day_sec, 2);
  aaa  = 1.;
  aaa2 = 1.;

  if(nx > 1) {
     for(int i = 0; i < nNode; i++)
         v0[0][i] /= aaa;
     for(int i = 0; i < nNode * nNode; i++) 
         Cv1v1[i] /= aaa2;
  }
  if(ny > 1) {
     for(int i = 0; i < nNode; i++)
         v0[1][i] /= aaa;
     for(int i = 0; i < nNode * nNode; i++)
         Cv2v2[i] /= aaa2;   
  }
  if(nz > 1) {
     for(int i = 0; i < nNode; i++)
         v0[2][i] /= aaa;
     for(int i = 0; i < nNode * nNode; i++)
         Cv3v3[i] /= aaa2;             
  }
  if(nx > 1 && ny > 1) {
     for(int i = 0; i < nNode * nNode; i++) 
         Cv1v2[i] /= aaa2;
  }
  if(nx > 1 && nz > 1) {
     for(int i = 0; i < nNode * nNode; i++) 
         Cv1v3[i] /= aaa2;
  }
  if(ny > 1 && nz > 1) {
     for(int i = 0; i < nNode * nNode; i++) 
         Cv2v3[i] /= aaa2;
  }
}

// volume balance and uncertainty balance at each cell.
void Velocity::check() {
  ofstream os1("vol.out", ios::out);
  os1 << setiosflags(ios::fixed | ios::showpoint) 
      << setprecision(5);
        
  if(nx > 1 && ny > 1) {
     int ijk, ijk_im1, ijk_jm1;
     int num_blocks = nNode;
     double vx1, vx2, vy1, vy2, poro, vol;
     double aaa  =   1.;
     if(unit == FIELD) aaa = stbPday_m3Psec;
     for(int k = 0; k < nz; ++k) {
         double dz2 = gridPtr->getBdz(k) * gridPtr->getBdz(k);
         for(int j = 0; j < ny; ++j) {
             double dy2        = gridPtr->getBdy(j) * gridPtr->getBdy(j);
             for(int i = 0; i < nx; ++i) {
                 double dx2 = gridPtr->getBdx(i) * gridPtr->getBdx(i);
                 double dxy = gridPtr->getBdx(i) * gridPtr->getBdy(j);          
                 ijk = gridPtr->getIndex( i, j, k );
                 poro = fluidPtr->getPoro(ijk);
                 vx2 = v0[0][ijk];
                 vy2 = v0[1][ijk];
                 if(i == 0) { 
                    vx1 = 0;
                 } else {
                    ijk_im1 = gridPtr->getIndex( i - 1, j, k );
                    vx1 = v0[0][ijk_im1];
                 }
                 if(j == 0) {
                    vy1 = 0;
                 }else {
                    ijk_jm1 = gridPtr->getIndex( i, j - 1, k );
                    vy1 = v0[1][ijk_jm1];
                 }
                 
                 vol = ( ( vx2 - vx1 ) * gridPtr->getBdy(j) 
                        +( vy2 - vy1 ) * gridPtr->getBdx(i)
                       ) * gridPtr->getBdz(k)  * poro / aaa;

                 // variance part
                 double Cvx2_vx2 = Cv1v1[ ijk + ijk * num_blocks ];
                 double Cvy2_vy2 = Cv2v2[ ijk + ijk * num_blocks ];
                 double Cvx2_vy2 = Cv1v2[ ijk + ijk * num_blocks ];
                 double Cvx1_vx1, Cvx1_vx2, Cvx2_vx1;
                 double Cvy1_vy1, Cvy1_vy2, Cvy2_vy1;
                 double Cvx1_vy1, Cvx1_vy2, Cvx2_vy1;
   
                 if(i == 0 && j == 0 ) {
                    Cvx1_vy1 = 0.;
                 } else {
                    ijk_im1 = gridPtr->getIndex( i - 1, j, k );
                    ijk_jm1 = gridPtr->getIndex( i, j - 1, k );
                    Cvx1_vy1 = Cv1v2[ ijk_im1 + ijk_jm1 * num_blocks ];
                 } 
                 if(i == 0 ) {
                    Cvx1_vx1 = 0.;
                    Cvx1_vx2 = 0.;
                    Cvx2_vx1 = 0.;
                    Cvx1_vy2 = 0.;      
                 } else {
                    ijk_im1 = gridPtr->getIndex( i - 1, j, k );         
                    Cvx1_vx1 = Cv1v1[ ijk_im1 + ijk_im1  * num_blocks ];
                    Cvx1_vx2 = Cv1v1[ ijk_im1 + ijk      * num_blocks ];
                    Cvx2_vx1 = Cv1v1[ ijk     + ijk_im1  * num_blocks ]; 
                    Cvx1_vy2 = Cv1v2[ ijk_im1 + ijk      * num_blocks ];      
                 }
                 if(j == 0) {
                    Cvy1_vy1 = 0.;
                    Cvy1_vy2 = 0.;
                    Cvy2_vy1 = 0.;
                    Cvx2_vy1 = 0.;
                 } else { 
                    ijk_jm1 = gridPtr->getIndex( i, j - 1, k );
                    Cvy1_vy1 = Cv2v2[ ijk_jm1 + ijk_jm1 * num_blocks ];
                    Cvy1_vy2 = Cv2v2[ ijk_jm1 + ijk     * num_blocks ];
                    Cvy2_vy1 = Cv2v2[ ijk     + ijk_jm1 * num_blocks ];
                    Cvx2_vy1 = Cv1v2[ ijk     + ijk_jm1 * num_blocks ];
                 }
                 double cqq_total =( ( Cvx2_vx2 + Cvx1_vx1 - Cvx1_vx2 - Cvx2_vx1 ) * dy2 
                                    +( Cvy2_vy2 + Cvy1_vy1 - Cvy1_vy2 - Cvy2_vy1 ) * dx2
                                    +( Cvx2_vy2 + Cvx1_vy1 - Cvx1_vy2 - Cvx2_vy1 ) * dxy * 2.
                                   ) * dz2 * poro * poro /aaa /aaa;
                 os1 << ijk <<' '
                     << vol <<' '
                     << cqq_total << endl;
                 /*
                 if(fabs(vol) > 1.0e-10 || fabs(cqq_total) > 1.0e-10 ) {
                    cout << ijk <<' '
                         << vol <<' '
                         << cqq_total << endl;
                 }*/
             }
         }
     }
  }
  os1.close();
  //exit(0); 
}


// volume balance and uncertainty balance at each cell.
void Velocity::check3d() {
  int ijk, ijk_im1, ijk_jm1, ijk_km1;
  int num_blocks = nNode;
  double vx1, vx2, vy1, vy2, vz1, vz2;
  double poro;
  double aaa  =   1.;
  if(unit == FIELD) aaa = stbPday_m3Psec;
        
  //1) Mean Part
  double* vol = new double [num_blocks];
  for(int i = 0; i < num_blocks; ++i) vol[i] = 0.; 

  double dx, dy, dz;
  if(nx > 1) {
     for(int k = 0; k < nz; ++k) {
         dz = gridPtr->getBdz(k);
         for(int j = 0; j < ny; ++j) {
             dy = gridPtr->getBdy(j);         
             for(int i = 0; i < nx; ++i) {
                 ijk = gridPtr->getIndex( i, j, k );
                 vx2 = v0[0][ijk];
                 if(i == 0) { 
                    vx1 = 0;
                 } else {
                    ijk_im1 = gridPtr->getIndex( i - 1, j, k );
                    vx1 = v0[0][ijk_im1];
                 }
                 vol[ijk] +=  ( vx2 - vx1 ) * dy * dz; 
             }
         }
     }
  }
  if(ny > 1) {
     for(int k = 0; k < nz; ++k) {
         dz = gridPtr->getBdz(k);             
         for(int j = 0; j < ny; ++j) {
             for(int i = 0; i < nx; ++i) {
                 dx = gridPtr->getBdx(i);
                 ijk = gridPtr->getIndex( i, j, k );
                 vy2 = v0[1][ijk];
                 if(j == 0) {
                    vy1 = 0;
                 }else {
                    ijk_jm1 = gridPtr->getIndex( i, j - 1, k );
                    vy1 = v0[1][ijk_jm1];
                 }
                 vol[ijk] += ( vy2 - vy1 ) * dx * dz;
             }
         }
     }
  }
  if(nz > 1) {
     for(int k = 0; k < nz; ++k) {
         for(int j = 0; j < ny; ++j) {
             dy = gridPtr->getBdy(j);
             for(int i = 0; i < nx; ++i) {
                 dx = gridPtr->getBdx(i);
                 ijk = gridPtr->getIndex( i, j, k );
                 vz2 = v0[2][ijk];
                 if(k == 0) {
                    vz1 = 0;
                 }else {
                    ijk_jm1 = gridPtr->getIndex( i, j, k - 1 );
                    vz1 = v0[2][ijk_jm1];
                 }
                 vol[ijk] += ( vz2 - vz1 ) * dx * dy;
             }
         }
     }
  }
  ofstream os1("vol_avg.out", ios::out);
  //os1 << setiosflags(ios::fixed | ios::showpoint) 
  //    << setprecision(5);
  for(int k = 0; k < nz; ++k) {
      for(int j = 0; j < ny; ++j) {
          for(int i = 0; i < nx; ++i) {
              ijk = gridPtr->getIndex( i, j, k );
              vol[ijk] *= fluidPtr->getPoro(ijk) / aaa;
              if(fabs(vol[ijk]) > 1.e-10)
                 os1 << i <<' '
                     << j <<' ' 
                     << k <<' '
                     << vol[ijk] << endl;
          }
      }
  }  
  os1.close();

  //2) Variance Part
  double dx2, dy2, dz2;
  double* vol_var = new double [num_blocks];
  for(int i = 0; i < num_blocks; ++i) vol_var[i] = 0.; 
  if(nx > 1) {
     double Cvx2_vx2, Cvx1_vx1, Cvx2_vx1;
     for(int k = 0; k < nz; ++k) {
         dz2 = gridPtr->getBdz(k) * gridPtr->getBdz(k);
         for(int j = 0; j < ny; ++j) {
             dy2 = gridPtr->getBdy(j) * gridPtr->getBdy(j);
             for(int i = 0; i < nx; ++i) {
                 ijk = gridPtr->getIndex( i, j, k );
                 Cvx2_vx2 = Cv1v1[ ijk + ijk * num_blocks ];
                 if(i == 0) { 
                    Cvx1_vx1 = Cvx2_vx1 = 0;
                 } else {
                    ijk_im1 = gridPtr->getIndex( i - 1, j, k );         
                    Cvx1_vx1 = Cv1v1[ ijk_im1 + ijk_im1  * num_blocks ];
                    Cvx2_vx1 = Cv1v1[ ijk     + ijk_im1  * num_blocks ]; 
                 }
                 vol_var[ijk] +=  ( Cvx2_vx2 + Cvx1_vx1 - 2. * Cvx2_vx1 ) * dy2 * dz2; 
             }
         }
     }
  }
  if(ny > 1) {
     double Cvy2_vy2, Cvy1_vy1, Cvy2_vy1;          
     for(int k = 0; k < nz; ++k) {
         dz2 = gridPtr->getBdz(k) * gridPtr->getBdz(k);
         for(int j = 0; j < ny; ++j) {
             for(int i = 0; i < nx; ++i) {
                 dx2 = gridPtr->getBdx(i) * gridPtr->getBdx(i);
                 ijk = gridPtr->getIndex( i, j, k );
                 Cvy2_vy2 = Cv2v2[ ijk + ijk * num_blocks ];
                 if(j == 0) {
                    Cvy1_vy1 = Cvy2_vy1 = 0.;
                 }else {
                    ijk_jm1 = gridPtr->getIndex( i, j - 1, k );
                    Cvy1_vy1 = Cv2v2[ ijk_jm1 + ijk_jm1 * num_blocks ];
                    Cvy2_vy1 = Cv2v2[ ijk     + ijk_jm1 * num_blocks ];
                 }
                 vol_var[ijk] += ( Cvy2_vy2 + Cvy1_vy1 - 2. * Cvy2_vy1 ) * dx2 * dz2;
             }
         }
     }
  }
  if(nz > 1) {
     double Cvz2_vz2, Cvz1_vz1, Cvz2_vz1;
     for(int k = 0; k < nz; ++k) {
         for(int j = 0; j < ny; ++j) {
             dy2 = gridPtr->getBdy(j) * gridPtr->getBdy(j);
             for(int i = 0; i < nx; ++i) {
                 dx2 = gridPtr->getBdx(i) * gridPtr->getBdx(i);
                 ijk = gridPtr->getIndex( i, j, k );
                 Cvz2_vz2 = Cv3v3[ ijk + ijk * num_blocks ];
                 if(k == 0) {
                    Cvz1_vz1 = Cvz2_vz1 = 0.;
                 }else {
                    ijk_km1 = gridPtr->getIndex( i, j, k - 1 );
                    Cvz1_vz1 = Cv3v3[ ijk_km1 + ijk_km1 * num_blocks ];
                    Cvz2_vz1 = Cv3v3[ ijk     + ijk_km1 * num_blocks ];
                 }
                 vol_var[ijk] += ( Cvz2_vz2 + Cvz1_vz1 - 2. * Cvz2_vz1 ) * dx2 * dy2;
             }
         }
     }
  }
  if(nx > 1 && ny > 1) {
     double Cvx2_vy2, Cvx1_vy1, Cvx1_vy2, Cvx2_vy1;
     for(int k = 0; k < nz; ++k) {
         dz2 = gridPtr->getBdz(k) * gridPtr->getBdz(k);
         for(int j = 0; j < ny; ++j) {
             dy = gridPtr->getBdy(j);         
             for(int i = 0; i < nx; ++i) {
                 dx = gridPtr->getBdx(i); 
                 ijk = gridPtr->getIndex( i, j, k );
                 Cvx2_vy2 = Cv1v2[ ijk + ijk * num_blocks ];
                 if(i == 0 || j == 0 ) {
                    Cvx1_vy1 = 0.;
                 } else {
                    ijk_im1 = gridPtr->getIndex( i - 1, j, k );
                    ijk_jm1 = gridPtr->getIndex( i, j - 1, k );
                    Cvx1_vy1 = Cv1v2[ ijk_im1 + ijk_jm1 * num_blocks ];
                 } 
                 if(i == 0 ) {
                    Cvx1_vy2 = 0.;      
                 } else {
                    ijk_im1 = gridPtr->getIndex( i - 1, j, k );         
                    Cvx1_vy2 = Cv1v2[ ijk_im1 + ijk      * num_blocks ];      
                 }
                 if(j == 0) {
                    Cvx2_vy1 = 0.;
                 } else { 
                    ijk_jm1 = gridPtr->getIndex( i, j - 1, k );
                    Cvx2_vy1 = Cv1v2[ ijk     + ijk_jm1 * num_blocks ];
                 }
                 vol_var[ijk] += 2. * ( Cvx2_vy2 + Cvx1_vy1 - Cvx2_vy1 - Cvx1_vy2) * dx * dy * dz2; 
             }
         }
     }
  }
  if(nx > 1 && nz > 1) {
     double Cvx2_vz2, Cvx1_vz1, Cvx1_vz2, Cvx2_vz1;
     for(int k = 0; k < nz; ++k) {
         dz = gridPtr->getBdz(k);
         for(int j = 0; j < ny; ++j) {
             dy2 = gridPtr->getBdy(j) * gridPtr->getBdy(j);
             for(int i = 0; i < nx; ++i) {
                 dx = gridPtr->getBdx(i); 
                 ijk = gridPtr->getIndex( i, j, k );
                 Cvx2_vz2 = Cv1v3[ ijk + ijk * num_blocks ];
                 if(i == 0 || k == 0 ) {
                    Cvx1_vz1 = 0.;
                 } else {
                    ijk_im1 = gridPtr->getIndex( i - 1, j, k );
                    ijk_km1 = gridPtr->getIndex( i, j, k  - 1);
                    Cvx1_vz1 = Cv1v3[ ijk_im1 + ijk_km1 * num_blocks ];
                 } 
                 if(i == 0 ) {
                    Cvx1_vz2 = 0.;      
                 } else {
                    ijk_im1 = gridPtr->getIndex( i - 1, j, k );         
                    Cvx1_vz2 = Cv1v3[ ijk_im1 + ijk      * num_blocks ];      
                 }
                 if(k == 0) {
                    Cvx2_vz1 = 0.;
                 } else { 
                    ijk_km1 = gridPtr->getIndex( i, j, k - 1);
                    Cvx2_vz1 = Cv1v3[ ijk     + ijk_km1 * num_blocks ];
                 }
                 vol_var[ijk] += 2. * ( Cvx2_vz2 + Cvx1_vz1 - Cvx2_vz1 - Cvx1_vz2) * dx * dy2 * dz; 
             }
         }
     }
  }
  if(ny > 1 && nz > 1) {
     double Cvy2_vz2, Cvy1_vz1, Cvy1_vz2, Cvy2_vz1;
     for(int k = 0; k < nz; ++k) {
         dz = gridPtr->getBdz(k);
         for(int j = 0; j < ny; ++j) {
             dy = gridPtr->getBdy(j);
             for(int i = 0; i < nx; ++i) {
                 dx2 = gridPtr->getBdx(i) * gridPtr->getBdx(i); 
                 ijk = gridPtr->getIndex( i, j, k );
                 Cvy2_vz2 = Cv2v3[ ijk + ijk * num_blocks ];
                 if(j == 0 || k == 0 ) {
                    Cvy1_vz1 = 0.;
                 } else {
                    ijk_jm1 = gridPtr->getIndex( i, j - 1, k );
                    ijk_km1 = gridPtr->getIndex( i, j, k  - 1);
                    Cvy1_vz1 = Cv2v3[ ijk_jm1 + ijk_km1 * num_blocks ];
                 } 
                 if(j == 0 ) {
                    Cvy1_vz2 = 0.;      
                 } else {
                    ijk_jm1 = gridPtr->getIndex( i, j - 1, k );         
                    Cvy1_vz2 = Cv2v3[ ijk_jm1 + ijk      * num_blocks ];      
                 }
                 if(k == 0) {
                    Cvy2_vz1 = 0.;
                 } else { 
                    ijk_km1 = gridPtr->getIndex( i, j, k - 1);
                    Cvy2_vz1 = Cv2v3[ ijk     + ijk_km1 * num_blocks ];
                 }
                 vol_var[ijk] += 2. * ( Cvy2_vz2 + Cvy1_vz1 - Cvy2_vz1 - Cvy1_vz2) * dx2 * dy * dz; 
             }
         }
     }
  }

  ofstream os2("vol_var.out", ios::out);
  //os2 << setiosflags(ios::fixed | ios::showpoint) 
  //    << setprecision(5);
  for(int k = 0; k < nz; ++k) {
      for(int j = 0; j < ny; ++j) {
          for(int i = 0; i < nx; ++i) {
              ijk = gridPtr->getIndex( i, j, k );
              vol_var[ijk] *= (   fluidPtr->getPoro(ijk) / aaa 
                                * fluidPtr->getPoro(ijk) / aaa
                              );
              if(fabs(vol_var[ijk]) > 1.e-10)
               os2 << i <<' '
                   << j <<' ' 
                   << k <<' '
                   << vol_var[ijk] << endl;
          }
      }
      os2 << endl;
  }  
  os2.close();
  delete [] vol_var;
  delete [] vol;
 
  return;
  //exit(0); 
  ofstream os3("vel_cov.out", ios::out);
  int ijk1, ijk2;
  for(int k2 = 0; k2 < nz; ++k2) {
      for(int j2 = 0; j2 < ny; ++j2) {
          for(int i2 = 0; i2 < nx; ++i2) {
              ijk2 = gridPtr->getIndex( i2, j2, k2);        
              for(int k1 = 0; k1 < nz; ++k1) {
                  for(int j1 = 0; j1 < ny; ++j1) {
                      for(int i1 = 0; i1 < nx; ++i1) {
                          ijk1 = gridPtr->getIndex( i1, j1, k1);
                          if(i1 == nx -1 || i2 == nx -1) {
                             if(fabs(Cv1v1[ijk1 +ijk2 * num_blocks]) > 1.e-10)
                                os3<<"Cv1v1 = "
                                   << Cv1v1[ijk1 +ijk2 * num_blocks]<<endl;
                          }
                          if(j1 == ny -1 || j2 == ny -1)
                             if(fabs(Cv2v2[ijk1 +ijk2 * num_blocks]) > 1.e-10)
                                os3<<"Cv2v2 = "
                                   << Cv2v2[ijk1 +ijk2 * num_blocks]<<endl;
                          if(k1 == nz -1 || k2 == nz -1)
                             if(fabs(Cv3v3[ijk1 +ijk2 * num_blocks]) > 1.e-10)
                                os3<<"Cv3v3 = "
                                   << Cv3v3[ijk1 +ijk2 * num_blocks]<<endl;          
                          if(i1 == nx -1 || j2 == ny -1)
                             if(fabs(Cv1v2[ijk1 +ijk2 * num_blocks]) > 1.e-10)
                                os3<<"Cv1v2 = "
                                   << Cv1v2[ijk1 +ijk2 * num_blocks]<<endl;          
                          if(i1 == nx -1 || k2 == nz -1)
                             if(fabs(Cv1v3[ijk1 +ijk2 * num_blocks]) > 1.e-10)
                                os3<<"Cv1v3 = "
                                   << Cv1v3[ijk1 +ijk2 * num_blocks]<<endl;
                          if(j1 == ny -1 || k2 == nz -1)
                             if(fabs(Cv2v3[ijk1 +ijk2 * num_blocks]) > 1.e-10)
                                os3<<"Cv2v3 = "
                                   << Cv2v3[ijk1 +ijk2 * num_blocks]<<endl;          
                      }
                  }
              } 
          }
      }
  } 
  
  os3.close();
    
  ofstream os4("vel_avg.out", ios::out);
  os4 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);

  for(int k2 = 0; k2 < nz; ++k2) {
      for(int j2 = 0; j2 < ny; ++j2) {
          for(int i2 = 0; i2 < nx; ++i2) {
              ijk2 = gridPtr->getIndex( i2, j2, k2);        
              os4<< i2 << ' '
                 << j2 << ' '
                 << k2 << ' '
                 <<v0[0][ijk2] <<' '
                 <<v0[1][ijk2] <<' '
                 <<v0[2][ijk2] <<endl;
          }
          os4<<endl;
      }
      os4<<endl;
  }

  ofstream os5("vel_cov.out", ios::out);
  for(int k2 = 0; k2 < nz; ++k2) {
      for(int j2 = 0; j2 < ny; ++j2) {
          for(int i2 = 0; i2 < nx; ++i2) {
              ijk2 = gridPtr->getIndex( i2, j2, k2);        
              for(int k1 = 0; k1 < nz; ++k1) {
                  for(int j1 = 0; j1 < ny; ++j1) {
                      for(int i1 = 0; i1 < nx; ++i1) {
                          ijk1 = gridPtr->getIndex( i1, j1, k1);
                           os5 << Cv1v1[ijk1 +ijk2 * num_blocks] << ' '
                               << Cv2v2[ijk1 +ijk2 * num_blocks] << ' '
                               << Cv3v3[ijk1 +ijk2 * num_blocks] <<endl;          
                      }
                  }
              } 
          }
      }
  } 
  
  os5.close();
   
  //exit(0);
}

//===== writeCv() =====
void Velocity::writeCv(const Control &c) {
  
  char file1[100], file2[100], file3[100], file4[100];
  int i_ref = c.getI_Ref();
  int j_ref = c.getJ_Ref();
  int k_ref = c.getK_Ref();
  cout << "i_ref = " << i_ref << ' '
       << "j_ref = " << j_ref << ' '
       << "k_ref = " << k_ref << ' '
       << endl;       

  double aaa  =     ft_m/day_sec ;
  double aaa2 = pow(ft_m/day_sec, 2);
  if(unit != FIELD) {
     aaa  = 1.;
     aaa2 = 1.;
  }
  
  int ijk_ref = gridPtr->getIndex( i_ref, j_ref, k_ref );;
  double vx_var_ref = Cv1v1[ijk_ref + ijk_ref * nNode];
  double vy_var_ref = Cv2v2[ijk_ref + ijk_ref * nNode];
  if(nx > 1 && ny > 1) {
     sprintf( file1, "%sVELO_XY_1D.out", directory);
     sprintf( file2, "%sVELO_XY_2D.out", directory);
     sprintf( file3, "%sVELO_COV_XDir_1D.out", directory);
     sprintf( file4, "%sVELO_COV_Diag_1D.out", directory);

     ofstream os1(file1, ios::out);
     ofstream os2(file2, ios::out);
     ofstream os3(file3, ios::out);
     ofstream os4(file4, ios::out);
     os1 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);
     os2 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);
     os3 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);
     os4 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);

     int k = 0;
     // 2D XY Plot
     for(int j = 0; j < ny; ++j) {
         for(int i = 0; i < nx; ++i) {
             int ijk  = gridPtr->getIndex( i, j, k );
             os2<< i<<' '<<j<<' ';
             os2<< v0[0][ijk]/aaa <<' '
                << v0[1][ijk]/aaa <<' ';
             os2<< Cv1v1[ijk + ijk * nNode]/aaa2 <<' '
                << Cv2v2[ijk + ijk * nNode]/aaa2 <<' '
                << Cv1v2[ijk + ijk * nNode]/aaa2 << endl;
         }
         os2<<endl;
     }
     os2.close();

     // 1D X-Dir Plot
     for(int i = 0; i < nx - 1; ++i) {
         int j = j_ref;
         int ijk  = gridPtr->getIndex( i, j, k );
         double vx_var = Cv1v1[ijk + ijk * nNode];
         double vy_var = Cv2v2[ijk + ijk * nNode];
         double tmp11  = sqrt( vx_var * vx_var_ref );
         double tmp22  = sqrt( vy_var * vy_var_ref );
         double tmp12  = sqrt( vx_var * vy_var_ref );
         double cor11  = Cv1v1[ijk + ijk_ref * nNode]/tmp11;
         double cor22  = Cv2v2[ijk + ijk_ref * nNode]/tmp22;
         double cor12  = Cv1v2[ijk + ijk_ref * nNode]/tmp12;
         if(tmp11 < 0.0000001) cor11 = 0;
         if(tmp22 < 0.0000001) cor22 = 0;
         if(tmp12 < 0.0000001) cor12 = 0;
         os3<< i <<' ';
         os3<< Cv1v1[ijk + ijk_ref * nNode]/aaa2 <<' '
            << Cv2v2[ijk + ijk_ref * nNode]/aaa2 <<' '
            << Cv1v2[ijk + ijk_ref * nNode]/aaa2 <<' '
            << cor11 <<' '
            << cor22 <<' '
            << cor12 <<' '
            << endl;
     }
    
     if(nx == ny) {
	// 1D - XY Diag Plot
        for(int i = 0; i < nx - 1; ++i) {
            int j = i;
            int ijk  = gridPtr->getIndex( i, j, k );
            os1<< i<<' ';
            os1<< v0[0][ijk]/aaa <<' '
               << v0[1][ijk]/aaa <<' ';
            os1<< Cv1v1[ijk + ijk * nNode]/aaa2 <<' '
               << Cv2v2[ijk + ijk * nNode]/aaa2 <<' '
               << Cv1v2[ijk + ijk * nNode]/aaa2 << endl;

            double vx_var = Cv1v1[ijk + ijk * nNode];
            double vy_var = Cv2v2[ijk + ijk * nNode];
            double tmp11  = sqrt( vx_var * vx_var_ref );
            double tmp22  = sqrt( vy_var * vy_var_ref );
            double tmp12  = sqrt( vx_var * vy_var_ref );
            double cor11  = Cv1v1[ijk + ijk_ref * nNode]/tmp11;
            double cor22  = Cv2v2[ijk + ijk_ref * nNode]/tmp22;
            double cor12  = Cv1v2[ijk + ijk_ref * nNode]/tmp12;
            if(tmp11 < 0.0000001) cor11 = 0;
            if(tmp22 < 0.0000001) cor22 = 0;
            if(tmp12 < 0.0000001) cor12 = 0;
            os4<< i <<' ';
            os4<< Cv1v1[ijk + ijk_ref * nNode]/aaa2 <<' '
               << Cv2v2[ijk + ijk_ref * nNode]/aaa2 <<' '
               << Cv1v2[ijk + ijk_ref * nNode]/aaa2 <<' '
               << cor11 <<' '
               << cor22 <<' '
               << cor12 <<' '
               << endl;
        }
     } else {
        int j = j_ref;             
        for(int i = 0; i < nx; ++i) {
            int ijk  = gridPtr->getIndex( i, j, k );
            os1<< i <<' ';
            os1<< v0[0][ijk]/aaa <<' '
               << v0[1][ijk]/aaa <<' ';
            os1<< Cv1v1[ijk + ijk * nNode]/aaa2 <<' '
               << Cv2v2[ijk + ijk * nNode]/aaa2 <<' '
               << Cv1v2[ijk + ijk * nNode]/aaa2 << endl;
            double vx_var = Cv1v1[ijk + ijk * nNode];
            double vy_var = Cv2v2[ijk + ijk * nNode];
            double tmp11  = sqrt( vx_var * vx_var_ref );
            double tmp22  = sqrt( vy_var * vy_var_ref );
            double tmp12  = sqrt( vx_var * vy_var_ref );
            double cor11  = Cv1v1[ijk + ijk_ref * nNode]/tmp11;
            double cor22  = Cv2v2[ijk + ijk_ref * nNode]/tmp22;
            double cor12  = Cv1v2[ijk + ijk_ref * nNode]/tmp12;
            if(tmp11 < 0.0000001) cor11 = 0;
            if(tmp22 < 0.0000001) cor22 = 0;
            if(tmp12 < 0.0000001) cor12 = 0;
            os3<< i <<' ';
            os3<< Cv1v1[ijk + ijk_ref * nNode]/aaa2 <<' '
               << Cv2v2[ijk + ijk_ref * nNode]/aaa2 <<' '
               << Cv1v2[ijk + ijk_ref * nNode]/aaa2 <<' '
               << cor11 <<' '
               << cor22 <<' '
               << cor12 <<' '
               << endl;


        }
     }
     os1.close();
     os3.close();
     os4.close();
  }
  
  if(nx > 1 && nz > 1) {
     ofstream os1("VELO_XZ_1D.out", ios::out);          
     ofstream os2("VELO_XZ_2D.out", ios::out);
     os1 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);
     os2 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);
     int j = 0;
     for(int k = 0; k < nz; ++k) {
         for(int i = 0; i < nx; ++i) {
             int ijk  = gridPtr->getIndex( i, j, k );
             os2<< i<<' '<<k<<' ';
             os2<< v0[0][ijk] <<' '
                << v0[2][ijk] <<' ';
             os2<< Cv1v1[ijk + ijk * nNode] <<' '
                << Cv3v3[ijk + ijk * nNode] <<' '
                << Cv1v3[ijk + ijk * nNode] << endl;
         }
         os2<<endl;
     }
     os2.close();
  
     if(nz == nx ) {
        for(int i = 0; i < nx; ++i) {
            int k = i;
            int ijk  = gridPtr->getIndex( i, j, k );
            os1<< i <<' ';
            os1<< v0[0][ijk] <<' '
               << v0[2][ijk] <<' ';
            os1<< Cv1v1[ijk + ijk * nNode] <<' '
               << Cv3v3[ijk + ijk * nNode] <<' '
               << Cv1v3[ijk + ijk * nNode] <<endl;
        }
     } else {
        int k = nz/2;
        for(int i = 0; i < nx; ++i) {
            int ijk  = gridPtr->getIndex( i, j, k );
            os1<< i <<' ';
            os1<< v0[0][ijk] <<' '
               << v0[2][ijk] <<' ';
            os1<< Cv1v1[ijk + ijk * nNode] <<' '
               << Cv3v3[ijk + ijk * nNode] <<' '
               << Cv1v3[ijk + ijk * nNode] <<endl;
        }
     }
     os1.close();
  }

  if(ny > 1 && nz > 1) {
     ofstream os1("VELO_YZ_1D.out", ios::out);          
     ofstream os2("VELO_YZ_2D.out", ios::out);
     os1 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);
     os2 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);
     int i = 0;
     for(int k = 0; k < nz; ++k) {
         for(int j = 0; j < ny; ++j) {
             int ijk  = gridPtr->getIndex( i, j, k );
             os2<<  j<<' '<<k<<' ';
             os2<< v0[1][ijk] <<' '
                << v0[2][ijk] <<' ';
             os2<< Cv2v2[ijk + ijk * nNode] <<' '
                << Cv3v3[ijk + ijk * nNode] <<' ' 
                << Cv2v3[ijk + ijk * nNode] << endl;
         }
         os2<<endl;
     }
     os2.close();

     if(ny == nz ) {
        for(int j = 0; j < ny; ++j) {
            int k = j;
            int ijk  = gridPtr->getIndex( i, j, k );
            os1<< j          <<' ';
            os1<< v0[1][ijk] <<' '
               << v0[2][ijk] <<' ';
            os1<< Cv2v2[ijk + ijk * nNode] <<' '
               << Cv3v3[ijk + ijk * nNode] <<' '
               << Cv2v3[ijk + ijk * nNode] << endl;
        }
     } else {
        int k = nz/2;
        for(int j = 0; j < ny; ++j) {
            int ijk  = gridPtr->getIndex( i, j, k );
            os1<< j         <<' ';
            os1<< v0[1][ijk] <<' '
               << v0[2][ijk] <<' ';
            os1<< Cv2v2[ijk + ijk * nNode] <<' '
               << Cv3v3[ijk + ijk * nNode] <<' '
               << Cv2v3[ijk + ijk * nNode] << endl;
        }
     }     
     os1.close();
  }

  if(nx > 1 && ny > 1 && nz > 1) {
     ofstream os2("VELO_Avg_3D.out", ios::out);
     os2 << setiosflags(ios::fixed | ios::showpoint) 
         << setprecision(5);
     for(int k = 0; k < nz; ++k) {
         for(int j = 0; j < ny; ++j) {
             for(int i = 0; i < nx; ++i) {         
                 int ijk  = gridPtr->getIndex( i, j, k );
                 os2<< i<< ' ' << j <<' '<< k <<' ';
                 os2<< v0[0][ijk] <<' '
                    << v0[1][ijk] <<' '
                    << v0[2][ijk] <<' '
		    << Cv1v1[ijk + ijk * nNode] <<' '
		    << Cv2v2[ijk + ijk * nNode] <<' '
		    << Cv3v3[ijk + ijk * nNode] << endl;
             }
             os2<<endl;
         }
     }
     os2.close();
  }
}

