/////////////////////////////////////////////////
//                                              //
//                 Liyong Li                    //
//      Reservoir Simulation Research Team      //
//        Chevron Petroleum Technology Co.      //
//          1300 Beach Blvd., Rm. 3166          //
//            La Habra, CA 90631-6374           //
//             Phone:(562) 694-7366             //
//              Fax:(562) 694-7565              //
//            Email liyl@chevron.com            //
//                                              //
//                 Version 1.                   //
//                Apr. 07, 1999                 //
//////////////////////////////////////////////////
//                                                //
// Domain.cpp: implementation of the Domain class. //
//                                                //


#include "Domain.h"

// Construction/Destruction
Domain::Domain(ifstream &is, char *Dir) {
   debug = false;
   //debug = true;
   readUnit( is );
   
   addGridObj( is , Dir);
   addBlocks();
   addRegnObj( is , unit );
   addFluidObj(is , unit );
   addWellArr( is );
   addContObj( is );
   calcVolume();
   printWellArc();
   
   //cout << "in debug ..." << endl;
   //exit(0);
}

/*
Domain::Domain(Grid* grid, int nwells, Well *wells) {
   debug = false;

   gridPtr = grid;
   initialize(gridPtr->getNx(), gridPtr->getNy() , nwells);
   buildBlocks(gridPtr);
   buildWellArcs(wells);
}
*/

Domain::~Domain(){
   if(debug) cout<<"Domain::~Domain()"<<endl;
   
   for(int i = 0; i < node_; ++i) 
       delete blocks[i];
   delete []blocks;
   
   for(int i = 0; i < nWells; ++i)
       delete wellArcArray[i];
   delete []wellArcArray;

   delete []Cv1v1;
   delete []Cv2v2;
   delete []Cv1v2;
   delete []dxCv1v1;
   delete []dyCv1v1;
   delete []dxCv2v2;
   delete []dyCv2v2;
   delete []dx1Cv1v2;
   delete []dx2Cv1v2;
   delete []dy1Cv1v2;
   delete []dy2Cv1v2;

   // delete Objects;
   delete [] wellArr;
   delete fluidPtr;
   delete contPtr;
   delete regnPtr;
   delete gridPtr;
}

void  Domain::calcVolume() {
  volume = 0;
   double *bdx, *bdy, *bdz;
   bdx = gridPtr->getBdx();
   bdy = gridPtr->getBdy();
   bdz = gridPtr->getBdz();
  for(int k = 0; k < nz_; k++)
      for(int j = 0; j < ny_; j++)
          for(int i = 0; i < nx_; i++) {
              int ijk = gridPtr->getIndex(i, j, k);
	      volume += bdx[i] * bdy[j] * bdz[k] * fluidPtr->getPoro(ijk);
	  }
}

void Domain::readUnit(ifstream &is) {
  // --- read the unit from input file --- 
  cout << endl << endl;
  cout << "===== UNIT Info =====" << endl;
  int i; 
  Junk(is); is >> i; 
  unit = (Unit) i;
  (unit == FIELD)? cout << "FIELD UNIT\n" : cout << "METRIC UNIT\n";
}

void Domain::addGridObj(ifstream &is, char *Dir) {
  // --- read and build Grid ---
  gridPtr = new Grid(is, unit, Dir);
  nx_ = gridPtr->getNx();
  ny_ = gridPtr->getNy();
  nz_ = gridPtr->getNz();
  node_ = gridPtr->getnNode();
  gridPtr->PrintGrid( cout );
}

void Domain::addRegnObj(ifstream &is, Unit unit) {
  regnPtr = new Region(is, gridPtr, unit);
}

void Domain::addFluidObj( ifstream &is, Unit unit ) {
  fluidPtr = new Fluid( is, gridPtr, unit);
}

void Domain::addWellArr(ifstream &is) {
  int i; 
  Junk(is); is >> nWells;
  wellArr = (nWells > 0) ? new Well[ nWells ] : NULL;
  for(i = 0; i < nWells; i++) wellArr[i].readData(is, unit);
  cout << endl;
  cout << "===== Wells Info =====" << endl;
  cout << " Number of Wells: " << nWells << endl;
  for(i = 0; i < nWells; i++) {
      cout << endl;
      cout << "  --- Well No. " << i + 1 << " ---\n" 
           << wellArr[i];
  }
  addWellArcs( );
}

void Domain::addContObj(ifstream &is) {
  contPtr = new Control (is, gridPtr, unit);
  bool timeDepend = ( fluidPtr->isCompressible() )? true : false; 
  if( !timeDepend ) contPtr->setSteady();
  nDt = contPtr->getNDt();
  setBoundBlk(); 
  cout << *contPtr;
}

void Domain::setBoundBlk() {
  int nx = gridPtr->getNx(); 
  int ny = gridPtr->getNy(); 
  int i_ind = 0; 
  int j_ind = 0; 
  bool isInjc; 
  bool isProd;
  
  for ( int i = 0; i < 6; i++ ){ 
     isInjc = contPtr->isInjcBound(i);
     isProd = contPtr->isProdBound(i); 
     if ( i == 0 ) {
         i_ind = 0; 
         for ( j_ind = 0; j_ind < ny; j_ind++ ) {
             boundBlk = getBlock( i_ind, j_ind );
             if (isInjc) {
                 boundBlk -> setInjcBound(); 
             }
             if (isProd) {
                 boundBlk -> setProdBound(); 
             }
         }
     }
     if ( i == 1 ) {
         i_ind = nx - 1; 
         for ( j_ind = 0; j_ind < ny; j_ind++ ) {
             boundBlk = getBlock( i_ind, j_ind ); 
             if (isInjc) {
                 boundBlk -> setInjcBound(); 
             }
             if (isProd) {
                 boundBlk -> setProdBound(); 
             }
         }
     }
     if ( i == 2 ) {
         j_ind = 0;
         for ( i_ind = 0; i_ind < nx; i_ind++ ) {
             boundBlk = getBlock( i_ind, j_ind ); 
             if (isInjc) {
                 boundBlk -> setInjcBound(); 
             } 
             if (isProd) { 
                 boundBlk -> setProdBound(); 
             } 
         } 
         
     }
     if ( i == 3 ) {
         j_ind = ny-1;
         for ( i_ind = 0; i_ind < nx; i_ind++ ) {
             boundBlk = getBlock( i_ind, j_ind ); 
             if (isInjc) {
                 boundBlk -> setInjcBound(); 
             } 
             if (isProd) { 
                 boundBlk -> setProdBound(); 
             } 
         } 
         
     }
  }
   
}

void Domain::addBlocks() {
   blocks   = new Block* [ node_ ];        
   double *x0, *y0, *z0, *bdx, *bdy, *bdz;
   bdx = gridPtr->getBdx();
   bdy = gridPtr->getBdy();
   bdz = gridPtr->getBdz();
   x0  = gridPtr->getX();
   y0  = gridPtr->getY();
   z0  = gridPtr->getZ();
   if(gridPtr->getGridType() == POINT_DIS) {
      double xx0, yy0;
      for(int k = 0; k < nz_; ++k) {
          for(int j = 0; j < ny_; ++j) {
              if      (j == 0      ) yy0 = y0[j];
              else if (j == ny_ - 1) yy0 = y0[j] - bdy[j];
              else                   yy0 = y0[j] - bdy[j]/2.;
              for(int i = 0; i < nx_; ++i) {
                  if     (i == 0      ) xx0 = x0[i];
                  else if(i == nx_ - 1) xx0 = x0[i] - bdx[i];
                  else                  xx0 = x0[i] - bdx[i]/2.;
                  Point pnt_lb( xx0  , yy0  , i, j);
                  Point pnt_ct( x0[i], y0[j], i, j);
                  blocks[ getIndex(i,j) ] = new Block(bdx[i], bdy[j], bdz[k], i, j, k, pnt_lb, pnt_ct);
              }
          }
      }
   } else {
      for(int k = 0; k < nz_; ++k) {
          for(int j = 0; j < ny_; ++j) {
              for(int i = 0; i < nx_; ++i) {
                  Point pnt_lb( x0[i] - bdx[i]/2., y0[j] - bdy[j]/2., i, j);
                  Point pnt_ct( x0[i]            , y0[j]            , i, j);
                  blocks[ getIndex(i,j) ] = new Block(bdx[i], bdy[j], bdz[k], i, j, k, pnt_lb, pnt_ct);
              }
          }
      }
   }
}
/*
void Domain::buildWellArcs(int i, int nstrl, bool isInjc, bool isProd, 
   double well_x, double well_y, double radius, double degree1, double degree2) {
   wellArcArray[i] = new WellArc(blocks, isInjc, isProd, nstrl, nx_, ny_, well_x, well_y, radius, 
                                 degree1, degree2);
}
*/
void Domain::addWellArcs() {
   wellArcArray = new WellArc*[ nWells ];
   for(int i = 0; i < nWells; ++i ) {
       for(int j = 0;   j < wellArr[i].getWellLen(); ++j) {
           int nstrl      = wellArr[i].getNumPoints();
           int well_i     = wellArr[i].getWellIJK(0,j);
           int well_j     = wellArr[i].getWellIJK(1,j);
           double degree1 = wellArr[i].getDegree(0);
           double degree2 = wellArr[i].getDegree(1);
           double radius  = wellArr[i].getRadius();
           double well_x  = wellArr[i].getWellXYZ(0,j);
           double well_y  = wellArr[i].getWellXYZ(1,j);
           bool isInjc    = wellArr[i].isInjc();
           bool isProd    = wellArr[i].isProd();
           if(debug) {
              cout << well_i << ' '<< well_j  <<' '<<nstrl <<' '<<nx_<<' '<<ny_<<endl;
              cout << well_x << ' '<< well_y  <<' '<<radius<<' '<<endl;
              cout << degree1<< ' '<< degree2 << endl;
              cout << isInjc << ' '<< isProd  << endl;
           }
	   wellBlk = getBlock(well_i, well_j);

           if(isProd) {
              //getBlock(well_i, well_j)->setProdWell();
	      wellBlk->setProdWell();
           }
           if(isInjc) {
              //getBlock(well_i, well_j)->setInjcWell();
	      wellBlk->setInjcWell();
           }

           /*wellArcArray[i] = new WellArc(blocks, isInjc, isProd, nstrl, nx_, ny_, i, well_i, well_j,
                                         well_x, well_y, radius, degree1, degree2);*/
	   wellArcArray[i] = new WellArc(blocks, wellBlk, isInjc, isProd, nstrl, nx_, ny_, i, well_i, well_j,
                                         well_x, well_y, radius, degree1, degree2);
	   
       }
   }
   //printWellArc();
   //cout << "in debugging" <<endl;
   //exit(0);
}

void Domain::setVelocity(double* vx, double* vy){
   double vx1, vx2, vy1, vy2;
   BType *bType = contPtr->getBType();
   double *bCond = contPtr->getBCond(); 
   for(int j = 0; j < ny_; ++j) {
       for(int i = 0; i < nx_; ++i) {

          if(i == 0) {
             if (bType[0] == CONST_RATE) { 
                 int ij = gridPtr->getIndex( i, j, 0 ); 
                 vx1 = bCond[0]/fluidPtr->getPoro(ij);
             } else {
                 vx1 = 0.0; 
             }
             vx2 = vx[i   + j * nx_];
          }
          else if ( i == nx_ -1 ) {
             vx1 = vx[i-1 + j * nx_];
             if (bType[1] == CONST_RATE) {
                 int ij = gridPtr->getIndex( i, j, 0 ); 
                 vx2 = bCond[1]/fluidPtr->getPoro(ij); 
             } else {
                 vx2 = 0.0; 
             }
          }
          else {
             vx1 = vx[i-1 + j * nx_];
             vx2 = vx[i   + j * nx_];
          }

          if(j == 0) {
             if (bType[2] == CONST_RATE) {
                 int ij = gridPtr->getIndex( i, j, 0 ); 
                 vy1 = bCond[2]/fluidPtr->getPoro(ij); 
             } else {
                 vy1 = 0.0;
             }
             vy2 = vy[i   + j * nx_];
          }
          else if ( j == ny_ - 1) {
             vy1 = vy[i   + (j - 1) * nx_];
             if (bType[3] == CONST_RATE) {
                 int ij = gridPtr->getIndex( i, j, 0 ); 
                 vy2 = bCond[3]/fluidPtr->getPoro(ij);  
             } else {
                 vy2 = 0.0; 
             }
          }
          else {
             vy1 = vy[i   + (j - 1) * nx_];
             vy2 = vy[i   +       j * nx_];
          }
          getBlock(i, j)->setVelocity(vx1, vx2, vy1, vy2);
       }
   }
}

/*
void Domain::buildBlocks(int i, int j, double x0, double y0, double dx, double dy,
   double vx1, double vx2, double vy1, double vy2) {
   Point pnt_lb( x0 - dx/2., y0 - dy/2., i, j);
   Point pnt_ct( x0, y0, i, j);
   blocks[ getIndex(i,j) ] = new Block(vx1, vx2, vy1, vy2, dx, dy, i, j, pnt_lb, pnt_ct);
}
*/

void Domain::addVeloCova(double* tmp11, double* tmp22, double* tmp12) {
   Cv1v1    = new double [ node_ * node_ ];
   Cv2v2    = new double [ node_ * node_ ];
   Cv1v2    = new double [ node_ * node_ ];
   for(int i = 0; i < node_ * node_; ++i) {
      Cv1v1[i] = tmp11[i];
      Cv2v2[i] = tmp22[i];
      Cv1v2[i] = tmp12[i];
   }
}

void Domain::deleteVeloCovariances()
{
	delete []Cv1v1;
	delete []Cv2v2;
	delete []Cv1v2;
	delete []dxCv1v1;
	delete []dyCv1v1;
	delete []dxCv2v2;
	delete []dyCv2v2;
	delete []dx1Cv1v2;
	delete []dx2Cv1v2;
	delete []dy1Cv1v2;
	delete []dy2Cv1v2;
}

void Domain::readMCVmoments(){
char* file1="/.automount/matrix/d/a/data/liyl/stochastic/monte_carlo/quarter2d/multi_region/test/Velo_Vari_2dxy_cell.out";
  ifstream in1;
  in1.open(file1);
  if( in1.bad() ) {
      cerr << " Can not open " <<  file1  << endl;
      exit(8);
  }
  double *vx_avg = new double[node_];
  double *vy_avg = new double[node_];
  for(int i = 0; i < node_; ++i) vx_avg[i] = vy_avg[i] = 0;

  int dummy1, dummy2;
  double dummy3,dummy4;
  for(int j = 0; j < ny_; ++j) {
     for(int i = 0; i < nx_; ++i) {
         in1 >> dummy1 >> dummy2>> vx_avg[i+j*nx_]>>vy_avg[i+j*nx_]
             >> dummy3 >> dummy4; 
         //cout<<vx_avg[i+j*nx_]<<' '<<vy_avg[i+j*nx_]<<endl;
     }
  }
  in1.close();
  setVelocity(vx_avg, vy_avg);
  calcDerivative();
  delete [] vx_avg;
  delete [] vy_avg;
  
char* file2="/.automount/matrix/d/a/data/liyl/stochastic/monte_carlo/quarter2d/multi_region/test/MCVeloCov.out";
  ifstream in2;
  in2.open(file2);
  if( in2.bad() ) {
      cerr << " Can not open " <<  file2  << endl;
      exit(8);
  }
  for(int j = 0; j < node_ ; ++j) {
      for(int i = 0; i < node_ ; ++i) {
          in2 >> Cv1v1[i + j * node_ ] >>
                 Cv2v2[i + j * node_ ] >>
                 Cv1v2[i + j * node_ ];
          /*
          cout<< Cv1v1[i + j * node_ ]<< ' '
              << Cv2v2[i + j * node_ ]<< ' '
              << Cv1v2[i + j * node_ ]<<endl;
              */
      }
  }
  in2.close();
}

void Domain::addWellObj(int well_i, int well_j, bool isProd) {
   if( isProd ) { 
      blocks[ getIndex(well_i, well_j)]->setWellBlock();
   }
}

Block * Domain::point2block(Point& pnt) {
   return blocks[ getIndex( pnt.getI(), pnt.getJ() ) ];
}

void Domain::changeBlockVeloSign(){
   for(int i = 0; i < node_; ++i) {
       blocks[i]->changeVeloSign();
   }
}

void Domain::changeBlockProdInjc(){
   for(int i = 0; i < node_; ++i) {
       blocks[i]->switchProdInjc();
   }
   for(int i = 0; i < nWells; ++i) {
       wellArcArray[i]->switchProdInjc();
   }
}

void Domain::PrintDomainInfo(ostream & os) {
   for(int i = 0; i < node_; ++i ) {
       blocks[i]->PrintBlock(os);
   }
}

void Domain::calcDerivative() {
   dxCv1v1  = new double [ node_ * node_ ];
   dyCv1v1  = new double [ node_ * node_ ];
   dxCv2v2  = new double [ node_ * node_ ];
   dyCv2v2  = new double [ node_ * node_ ];
   dx1Cv1v2 = new double [ node_ * node_ ];
   dx2Cv1v2 = new double [ node_ * node_ ];
   dy1Cv1v2 = new double [ node_ * node_ ];
   dy2Cv1v2 = new double [ node_ * node_ ];
        
   int node1, node2;
   int node2_xm1, node2_ym1, node1_xm1, node1_ym1;
   double *bdx, *bdy;
   bdx = gridPtr->getBdx();
   bdy = gridPtr->getBdy();
   for(int j2 = 0; j2 < ny_; ++j2){
       for(int i2 = 0; i2 < nx_; ++i2) {
           node2     = i2     +  j2      * nx_;
           node2_xm1 = i2 - 1 +  j2      * nx_;
           node2_ym1 = i2     + (j2 - 1) * nx_;
           for(int j1 = 0; j1 < ny_; ++j1){        
               for(int i1 = 0; i1 < nx_; ++i1) {
                   node1     = i1     +  j1      * nx_;
                   node1_xm1 = i1 - 1 +  j1      * nx_;
                   node1_ym1 = i1     + (j1 - 1) * nx_;
                   dx2Cv1v2[node1 + node2 * node_] = 0.;
                   dyCv1v1 [node1 + node2 * node_] = 0.;
                   dxCv2v2 [node1 + node2 * node_] = 0.;
                   dy1Cv1v2[node1 + node2 * node_] = 0.;
                   if( i1 == 0) {
                     dxCv1v1 [node1 + node2 * node_] =
                                Cv1v1[node1     + node2 * node_]  /bdx[i1];
                     dx1Cv1v2[node1 + node2 * node_] =
                                Cv1v2[node1     + node2 * node_]  /bdx[i1];
                   } else {
                     dxCv1v1 [node1 + node2 * node_] =
                             (  Cv1v1[node1     + node2 * node_]
                              - Cv1v1[node1_xm1 + node2 * node_] )/bdx[i1];
                     dx1Cv1v2[node1 + node2 * node_] =
                             (  Cv1v2[node1     + node2 * node_]
                              - Cv1v2[node1_xm1 + node2 * node_] )/bdx[i1];
                   }           
                   if( j1 == 0) {
                     dyCv2v2 [node1 + node2 * node_] = 
                                Cv2v2[node1     + node2 * node_]  /bdy[j1];
                     dy2Cv1v2[node2 + node1 * node_] = 
                                Cv1v2[node2 + node1     * node_]  /bdy[j1];
                   } else {
                     dyCv2v2 [node1 + node2 * node_] = 
                             (  Cv2v2[node1     + node2 * node_] 
                              - Cv2v2[node1_ym1 + node2 * node_] )/bdy[j1];
                     dy2Cv1v2[node2 + node1 * node_] = 
                             (  Cv1v2[node2 + node1     * node_] 
                              - Cv1v2[node2 + node1_ym1 * node_] )/bdy[j1];
                   }
               }
           }
       }
   }
}

//
void Domain::printWellArc() {
   ofstream fileOu ("WellArc.out", ios::out);
   for(int i = 0; i < nWells; ++i ) {
       wellArcArray[i]->PrintArc(fileOu);
   }
   fileOu.close();
}
