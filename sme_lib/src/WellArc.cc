#include "WellArc.h"

WellArc::WellArc(Block** blocks, Block* wellBlock, bool injc, bool prod, int len, int nx, int ny, int iw, int i0, int j0, 
                 double x0, double y0, double radius, double angle1 = 0., double angle2 = Pi/2.)
:blockPtr(blocks), wellBlockPtr(wellBlock), wellType_Injc(injc), wellType_Prod(prod), len_(len),  radius_(radius),
id(iw), i_(i0), j_(j0) {
   cout << "WellArc::WellArc()"<< endl;
   angle1 *= Pi/180.;
   angle2 *= Pi/180.;

   smallValue = 1.0e-8;
   // allocate memory
   allocMem (); 
   
   // Injector only!
   if(wellType_Injc == 1) { 

      int flag_scheme = 3;
   
      if( flag_scheme == 1) {
          double total_degree = angle2 - angle1;
          double dtheta = total_degree / (len_);
          for(int i = 0; i < len_; i++){
              theta[i]  = angle1 + dtheta/2. + dtheta * i;
              weight[i] = dtheta;
              theta[i]  *= radius;     
          }
          CalcCirclePoint(nx, ny, x0, y0);
      } else if(flag_scheme == 2) {
          gauleg(angle1, angle2, theta, weight, len_);
          CalcCirclePoint(nx, ny, x0, y0);
          for(int i = 0; i < len_; i++){
              theta[i]  *= radius;     
          }
      } else if(flag_scheme == 3){
          FaceType(nx, ny, i0, j0);
      } else {
         cout<<"not implemented yet! " << endl;
         exit(0);
      }
   }
   
}

void WellArc::allocMem () {
   ptSet_ = new Point [ len_ ];
   theta  = new double[ len_ ];
   weight = new double[ len_ ];
   v_total= new double[ len_ ];
   //vx     = new double[ len_ ];
   //vy     = new double[ len_ ];
}

void WellArc::free() { 
  delete [] ptSet_; 
  delete [] theta;
  delete [] weight;
  delete [] v_total;
  //delete [] vx;
  //delete [] vy;
}

WellArc& WellArc::operator=(const WellArc & a) {
  if (this != &a) {
      free();
      copy(a);
  } 
  return *this;
}

void WellArc::copy(const WellArc& rhs) {
     ptSet_ = new Point[len_ = rhs.len_];
     for (int i = 0; i < len_; i++)
          ptSet_[i] = rhs.ptSet_[i];
}

Point& WellArc::operator[](int i) {
     return ptSet_[i];
}

const Point& WellArc::operator[](int i) const {
     return ptSet_[i];
}

void WellArc::switchProdInjc(){
     bool tmp = wellType_Injc;
     wellType_Injc = wellType_Prod;
     wellType_Prod = tmp;
}

void WellArc::setPoint(int i, int iblock, int jblock, double xtemp, double ytemp) {
   ptSet_[i] = Point (xtemp, ytemp, iblock, jblock );
}

//arrange Point Set Order 
void WellArc::setPointOrder(int nx, int ny) {
   //cout << "WellArc::setPointOrder(Block* thiswellBlockPtr, int nx, int ny)" << endl;

   double cell_xlen = wellBlockPtr->getDx(); //thiswellBlockPtr->getDx();
   double cell_ylen = wellBlockPtr->getDy(); //thiswellBlockPtr->getDy();
   double cell_zlen = wellBlockPtr->getDz(); //()thiswellBlockPtr->getDz();
   double x1        = wellBlockPtr->getBlockLbPnt().getX(); // thiswellBlockPtr->getBlockLbPnt().getX();
   double y1        = wellBlockPtr->getBlockLbPnt().getY(); //thiswellBlockPtr->getBlockLbPnt().getY();
   double x2 = x1 + cell_xlen;
   double y2 = y1 + cell_ylen; 

  // cout << "x: " << x1 <<' ' << x2 << endl;
  // cout << "y: " << y1 <<' ' << y2 << endl;

   int *num_face   = new int[4];
   int *index_face = new int[len_];

   for(int i = 0; i < 4; ++i) num_face[i] = 0;

   for(int i = 0; i < len_ ; ++i) {
       if( fabs(ptSet_[i].getX() - x1) < smallValue) {++num_face[0]; index_face[i] = 0; }
       if( fabs(ptSet_[i].getY() - y1) < smallValue) {++num_face[1]; index_face[i] = 1; }
       if( fabs(ptSet_[i].getX() - x2) < smallValue) {++num_face[2]; index_face[i] = 2; }
       if( fabs(ptSet_[i].getY() - y2) < smallValue) {++num_face[3]; index_face[i] = 3; }
   } 
   
   for(int i = 0; i < len_; ++i) weight[i] = 1.;

   double *dx1, *dx2, *dy1, *dy2;
   dx1 = (num_face[0] > 0) ? new double [ num_face[0] + 1 ] : NULL;
   dy1 = (num_face[1] > 0) ? new double [ num_face[1] + 1 ] : NULL;
   dx2 = (num_face[2] > 0) ? new double [ num_face[2] + 1 ] : NULL;
   dy2 = (num_face[3] > 0) ? new double [ num_face[3] + 1 ] : NULL;
   
   double tmp, ds1, ds2;
   int im1, ip1;
   for(int i = 0; i < len_; ++i) {
       im1 = i - 1;
       ip1 = i + 1;

       if( i == 0       ) im1 = len_ - 1;
       if( i == len_ - 1) ip1 = 0;
      
       // along x1 or x2 lines
       if(index_face[i] == 0 || index_face[i] == 2) {
          if(index_face[im1] == index_face[i] ) {
             ds1 = ptSet_[im1].getY();
          } else {
             ds1 = y1;
             if( fabs (y2 - ptSet_[i].getY() ) < fabs(y1 - ptSet_[i].getY()) ) 
             ds1 = y2;
          }
          if(index_face[ip1] == index_face[i] ) {
             ds2 = ptSet_[ip1].getY();
          } else {
             ds2 = y1;
             if( fabs (y2 - ptSet_[i].getY() ) < fabs(y1 - ptSet_[i].getY())  ) 
             ds2 = y2;
          }

          if(index_face[im1] != index_face[i] ) {
             theta[i] = fabs(ptSet_[i].getY() - ds1)
                      + fabs(ds2 - ptSet_[i].getY() ) /2.;
          } else if(index_face[ip1] != index_face[i] ) {
             theta[i] = fabs(ptSet_[i].getY() - ds1) /2.
                      + fabs(ds2 - ptSet_[i].getY() );
          } else {
            theta[i] = fabs(ds2 - ds1 )/2.;
          }
       }

       if(index_face[i] == 1 || index_face[i] == 3) {
          if(index_face[im1] == index_face[i] ) {
             ds1 = ptSet_[im1].getX();
          } else {
             ds1 = x1;
             if( fabs (x2 - ptSet_[i].getX() ) < fabs(x1 - ptSet_[i].getX())  ) 
             ds1 = x2;
          }
          if(index_face[ip1] == index_face[i] ) {
             ds2 = ptSet_[ip1].getX();
          } else {
             ds2 = x1;
             if( fabs (x2 - ptSet_[i].getX() ) < fabs(x1 - ptSet_[i].getX())  ) 
             ds2 = x2;
          }
          if(index_face[im1] != index_face[i] ) {
             theta[i] = fabs(ptSet_[i].getX() - ds1)
                      + fabs(ds2 - ptSet_[i].getX() ) /2.;
          } else if(index_face[ip1] != index_face[i] ) {
             theta[i] = fabs(ptSet_[i].getX() - ds1) /2.
                      + fabs(ds2 - ptSet_[i].getX() );
          } else {
            theta[i] = fabs(ds2 - ds1 )/2.;
          }
       }
   }
   
   double sum0, sum1, sum2, sum3;
   sum0 = sum1 = sum2 = sum3 = 0.;
   for(int i = 0; i < len_ ; ++i) {
       if(index_face[i] == 0) sum0 += theta[i];
       if(index_face[i] == 1) sum1 += theta[i];  
       if(index_face[i] == 2) sum2 += theta[i];
       if(index_face[i] == 3) sum3 += theta[i];
       /*
       cout << index_face[i] << ' '
            << ptSet_[i].getX() << ' '
            << ptSet_[i].getY() << ' '
            << theta[i] << ' '
            << weight[i] << ' '
            << endl; */
       
   } 
       //thiswellBlockPtr
   /*
   cout << "sum0 = " <<  sum0 << endl;
   cout << "sum1 = " <<  sum1 << endl;
   cout << "sum2 = " <<  sum2 << endl;
   cout << "sum3 = " <<  sum3 << endl;
   */
   if( wellBlockPtr->getI() != 0 && fabs(sum0 - cell_ylen) > 0.0000001 ) {
       cout << endl <<
       cout << "The calculation of theta might be wrong! " << endl;
       cout << "Check it out!" << endl;
       cout << "sum0 = " <<  sum0 << endl;
       exit(0);
   }
   if( wellBlockPtr->getJ() != 0 && fabs(sum1 - cell_xlen) > 0.0000001 ) {
       cout << endl <<
       cout << "The calculation of theta might be wrong! " << endl;
       cout << "Check it out!" << endl;
   
       cout << "sum1 = " <<  sum1 << endl;
       exit(0);
   }
   if( wellBlockPtr->getI() != nx - 1 && fabs(sum2 - cell_ylen) > 0.0000001 ) {
       cout << endl <<
       cout << "The calculation of theta might be wrong! " << endl;
       cout << "Check it out!" << endl;
       cout << "sum2 = " <<  sum2 << endl;
       exit(0);
   }
   if( wellBlockPtr->getJ() != ny - 1 && fabs(sum3 - cell_xlen) > 0.0000001 ) {
       cout << endl <<
       cout << "The calculation of theta might be wrong! " << endl;
       cout << "Check it out!" << endl;
       cout << "sum3 = " <<  sum3 << endl;
       exit(0);
   }
   //cout << sum0 <<' '<< sum1 <<' ' << sum2 << ' ' << sum3 << endl;
   //cout << endl << "!!!! OK !!!! " << endl;
   double cell_z_len = wellBlockPtr->getDz();
   for(int i = 0; i < len_; ++i) {
       if(index_face[i] == 0 ) {
          v_total[i] = wellBlockPtr->getVx1();
          theta[i]   =   theta[i] * cell_z_len;
       }
       if(index_face[i] == 1 ) {
          v_total[i] = wellBlockPtr->getVy1();
          theta[i]   =   theta[i] * cell_z_len;
       }
       if(index_face[i] == 2 ) {
          v_total[i] = wellBlockPtr->getVx2();
	  theta[i]   = -theta[i] * cell_z_len;
       }
       if(index_face[i] == 3 ) {
          v_total[i] = wellBlockPtr->getVy2();
	  theta[i]   = -theta[i] * cell_z_len;
       }
       //cout  << "vtotal = " << v_total[i] << endl;
   }
   //exit(0);
   double sum = 0.;
   for(int i = 0; i < len_; ++i) { 
       //cout << index_face[i]<<' '
       //     << theta[i] <<' '
       //     << v_total[i] << ' '
       //     << endl;    
       sum += v_total[i] * theta[i];
   }
   cout << "vtotal = " << sum <<' '<< sum / stbPday_m3Psec << endl;
   //exit(0);
   
   delete [] dx1;
   delete [] dx2;
   delete [] dy1;
   delete [] dy2;
   delete [] num_face;
   delete [] index_face;
}

// Streamline Launch Type 1
void WellArc::CalcCirclePoint(int nx, int ny, double x0, double y0) {
   double xtemp, ytemp;
   double x1, x2, y1, y2;
   int iblock, jblock;
      
   double sum = 0.;
   for(int i = 0; i < len_; i++){
       xtemp = radius_ * cos(theta[i]) + x0;
       ytemp = radius_ * sin(theta[i]) + y0;
       for(int itmp = 0; itmp < nx; ++itmp) {
           for(int jtmp = 0; jtmp < ny; ++jtmp) {
               x1 = blockPtr[itmp + jtmp * nx]->getBlockLbPnt().getX();
               y1 = blockPtr[itmp + jtmp * nx]->getBlockLbPnt().getY();               
               x2 = x1 + blockPtr[itmp + jtmp * nx]->getDx();
               y2 = y1 + blockPtr[itmp + jtmp * nx]->getDy();
               if(xtemp <= x2 && xtemp >= x1 && ytemp <= y2 && ytemp >= y1) {
                  iblock = itmp;
                  jblock = jtmp;
               }
           }
       }
       Point tmp_pt(xtemp, ytemp, iblock, jblock );
       ptSet_[i] = tmp_pt;
   }
}

void WellArc::FaceType(int nx, int ny, int i0, int j0){
   int nfacex1, nfacey1, nfacex2, nfacey2;
   
   nfacex1 = nfacey1 = nfacex2 = nfacey2 = 1;

   if(i0 == 0     ) nfacex1 = 0;
   if(i0 == nx - 1) nfacex2 = 0;
   if(j0 == 0     ) nfacey1 = 0;
   if(j0 == ny - 1) nfacey2 = 0;

   int nface = nfacex1 + nfacey1 + nfacex2 + nfacey2 ;
   
   int *num_face = new int[nface];
   num_face[0] = len_ / nface;
   for(int i = 1; i < nface; ++i) {
       int sum = 0;
       for(int j = 0; j < i; ++j) {
           sum += num_face[j];
       }
       num_face[i] = (len_ - sum)/(nface - i);
   }
   
   int sum = 0;
   for(int i = 0; i < nface; ++i) {
       sum += num_face[i];
       cout << "num_face = " << num_face[i] << endl;
   }
   if(sum != len_) {
      for(int i = 0; i < nface; ++i) 
         cout << num_face[i] << endl;
      cout << "the number of particle launching points is wrong!" << endl;
      exit(0);
   }
   
   double cell_xlen = wellBlockPtr->getDx();
   double cell_ylen = wellBlockPtr->getDy();
   double x1 = wellBlockPtr->getBlockLbPnt().getX();
   double y1 = wellBlockPtr->getBlockLbPnt().getY();
   double x2 = x1 + cell_xlen;
   double y2 = y1 + cell_ylen;

   int  count = -1;
   int istart =  0;
   int iend   =  0;

   if(nfacex2 == 1) {
         ++count;
         iend += num_face[count];
         double delta_y = cell_ylen / num_face[count];
         double ytemp = y2 - delta_y / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(x2, ytemp - delta_y * (i - istart), i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }     

      if(nfacey1 == 1) {
         ++count;
         iend += num_face[count];
         double delta_x = cell_xlen / num_face[count];
         double xtemp = x2 - delta_x / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(xtemp - delta_x * (i - istart), y1, i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }    
      if(nfacex1 == 1) {
         ++count;
         iend += num_face[count];
         double delta_y = cell_ylen / num_face[count];
         double ytemp = y1 + delta_y / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(x1, ytemp + delta_y * (i - istart), i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }
      if(nfacey2 == 1) {
         ++count;
         iend += num_face[count];
         double delta_x = cell_xlen / num_face[count];
         double xtemp = x1 + delta_x / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(xtemp + delta_x * (i - istart), y2, i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }

  /* 
   if(nfacex2 == 0 || nfacey2 == 0 ) {
      if(nfacex2 == 1) {
         ++count;
         iend += num_face[count];
         double delta_y = cell_ylen / num_face[count];
         double ytemp = y2 - delta_y / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(x2, ytemp - delta_y * (i - istart), i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }     

      if(nfacey1 == 1) {
         ++count;
         iend += num_face[count];
         double delta_x = cell_xlen / num_face[count];
         double xtemp = x2 - delta_x / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(xtemp - delta_x * (i - istart), y1, i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }    
      if(nfacex1 == 1) {
         ++count;
         iend += num_face[count];
         double delta_y = cell_ylen / num_face[count];
         double ytemp = y1 + delta_y / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(x1, ytemp + delta_y * (i - istart), i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }
      if(nfacey2 == 1) {
         ++count;
         iend += num_face[count];
         double delta_x = cell_xlen / num_face[count];
         double xtemp = x1 + delta_x / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(xtemp + delta_x * (i - istart), y2, i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }

   } else {
      if(nfacex1 == 1) {
         ++count;
         iend += num_face[count];
         double delta_y = cell_ylen / num_face[count];
         double ytemp = y1 + delta_y / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(x1, ytemp + delta_y * (i - istart), i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }
      if(nfacey2 == 1) {
         ++count;
         iend += num_face[count];
         double delta_x = cell_xlen / num_face[count];
         double xtemp = x1 + delta_x / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(xtemp + delta_x * (i - istart), y2, i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }
      if(nfacex2 == 1) {
         ++count;
         iend += num_face[count];
         double delta_y = cell_ylen / num_face[count];
         double ytemp = y2 - delta_y / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(x2, ytemp - delta_y * (i - istart), i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }     
      if(nfacey1 == 1) {
         ++count;
         iend += num_face[count];
         double delta_x = cell_xlen / num_face[count];
         double xtemp = x2 - delta_x / 2.;
         for(int i = istart; i < iend; ++i) {
             Point tmp_pt(xtemp - delta_x * (i - istart), y1, i0, j0);
             ptSet_[i] = tmp_pt;
         }
         istart += num_face[count];
      }    
   }
*/
   delete[] num_face;
}

// Print Out
void WellArc::PrintArc(ostream & os) {
   for(int i = 0; i < len_; i++)
       ptSet_[i].PrintXY(os);
/*
   ofstream fileOu ("CellBoundary.out", ios::out);
   for(int i = 0; i < len_; ++i) {
       fileOu << ptSet_[i].getX() <<' '
              << ptSet_[i].getY() <<endl;
   }
   fileOu.close();
   exit(0);
*/
}

void WellArc::DisplayArc() {
   for(int i = 0; i < len_; i++) {
       cout<< theta [i] << ' '
	   << weight[i] << ' '    
           << v_total[i];
       ptSet_[i].PrintXY(cout);
   }
}

